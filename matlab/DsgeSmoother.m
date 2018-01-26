function [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T,R,P,PK,decomp,trend_addition,state_uncertainty,M_,oo_,options_,bayestopt_] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_)
% Estimation of the smoothed variables and innovations.
%
% INPUTS
%   o xparam1       [double]   (p*1) vector of (estimated) parameters.
%   o gend          [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o data          [double]   (n*T) matrix of data.
%   o data_index    [cell]      1*smpl cell of column vectors of indices.
%   o missing_value 1 if missing values, 0 otherwise
%   o M_            [structure] decribing the model
%   o oo_           [structure] storing the results
%   o options_      [structure] describing the options
%   o bayestopt_    [structure] describing the priors
%   o estim_params_ [structure] characterizing parameters to be estimated
%
% OUTPUTS
%   o alphahat      [double]  (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
%   o etahat        [double]  (r*T) matrix, smoothed structural shocks (r>=n is the number of shocks).
%   o epsilonhat    [double]  (n*T) matrix, smoothed measurement errors.
%   o ahat          [double]  (m*T) matrix, updated (endogenous) variables (a_{t|t}) (decision-rule order)
%   o SteadyState   [double]  (m*1) vector specifying the steady state level of each endogenous variable (declaration order)
%   o trend_coeff   [double]  (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
%   o aK            [double]  (K,n,T+K) array, k (k=1,...,K) steps ahead
%                                   filtered (endogenous) variables  (decision-rule order)
%   o T and R       [double]  Matrices defining the state equation (T is the (m*m) transition matrix).
%   o P:            (m*m*(T+1)) 3D array of one-step ahead forecast error variance
%                       matrices (decision-rule order)
%   o PK:           (K*m*m*(T+K)) 4D array of k-step ahead forecast error variance
%                       matrices (meaningless for periods 1:d) (decision-rule order)
%   o decomp        (K*m*r*(T+K)) 4D array of shock decomposition of k-step ahead
%                       filtered variables (decision-rule order)
%   o trend_addition [double] (n*T) pure trend component; stored in options_.varobs order
%   o state_uncertainty [double] (K,K,T) array, storing the uncertainty
%                                   about the smoothed state (decision-rule order)
%   o M_            [structure] decribing the model
%   o oo_           [structure] storing the results
%   o options_      [structure] describing the options
%   o bayestopt_    [structure] describing the priors
%
% Notes:
%   m:  number of endogenous variables (M_.endo_nbr)
%   T:  number of Time periods (options_.nobs)
%   r:  number of strucural shocks (M_.exo_nbr)
%   n:  number of observables (length(options_.varobs))
%   K:  maximum forecast horizon (max(options_.nk))
%
%   To get variables that are stored in decision rule order in order of declaration
%   as in M_.endo_names, ones needs code along the lines of:
%   variables_declaration_order(dr.order_var,:) = alphahat
%
%   Defines bayestopt_.mf = bayestopt_.smoother_mf (positions of observed variables
%   and requested smoothed variables in decision rules (decision rule order)) and
%   passes it back via global variable
%
% ALGORITHM
%   Diffuse Kalman filter (Durbin and Koopman)
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2006-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

alphahat        = [];
etahat  = [];
epsilonhat      = [];
ahat          = [];
SteadyState   = [];
trend_coeff   = [];
aK            = [];
T             = [];
R             = [];
P             = [];
PK            = [];
decomp        = [];
vobs            = length(options_.varobs);
smpl          = size(Y,2);

if ~isempty(xparam1) %not calibrated model
    M_ = set_all_parameters(xparam1,estim_params_,M_);
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------
oldoo.restrict_var_list = oo_.dr.restrict_var_list;
oldoo.restrict_columns = oo_.dr.restrict_columns;
oo_.dr.restrict_var_list = bayestopt_.smoother_var_list;
oo_.dr.restrict_columns = bayestopt_.smoother_restrict_columns;

[T,R,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);

if info~=0
    print_info(info,options_.noprint, options_);
    return
end
oo_.dr.restrict_var_list = oldoo.restrict_var_list;
oo_.dr.restrict_columns = oldoo.restrict_columns;

%get location of observed variables and requested smoothed variables in
%decision rules
bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
if options_.noconstant
    constant = zeros(vobs,1);
else
    if options_.loglinear
        constant = log(SteadyState(bayestopt_.mfys));
    else
        constant = SteadyState(bayestopt_.mfys);
    end
end
trend_coeff = zeros(vobs,1);
if bayestopt_.with_trend == 1
    [trend_addition, trend_coeff] =compute_trend_coefficients(M_,options_,vobs,gend);
    trend = constant*ones(1,gend)+trend_addition;
else
    trend_addition=zeros(size(constant,1),gend);
    trend = constant*ones(1,gend);
end
start = options_.presample+1;
np    = size(T,1);
mf    = bayestopt_.mf;
% ------------------------------------------------------------------------------
%  3. Initial condition of the Kalman filter
% ------------------------------------------------------------------------------
%
%  Here, Pinf and Pstar are determined. If the model is stationary, determine
%  Pstar as the solution of the Lyapunov equation and set Pinf=[] (Notation follows
%  Koopman/Durbin (2003), Journal of Time Series Analysis 24(1))
%
Q = M_.Sigma_e;
H = M_.H;

if isequal(H,0)
    H = zeros(vobs,vobs);
end

Z = zeros(vobs,size(T,2));
for i=1:vobs
    Z(i,mf(i)) = 1;
end

expanded_state_vector_for_univariate_filter=0;
kalman_algo = options_.kalman_algo;
if options_.lik_init == 1               % Kalman filter
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    Pstar=lyapunov_solver(T,R,Q,options_);
    Pinf        = [];
elseif options_.lik_init == 2           % Old Diffuse Kalman filter
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    Pstar = options_.Harvey_scale_factor*eye(np);
    Pinf        = [];
elseif options_.lik_init == 3           % Diffuse Kalman filter
    if kalman_algo ~= 4
        kalman_algo = 3;
    else
        if ~all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                                                %Augment state vector (follows Section 6.4.3 of DK (2012))
            expanded_state_vector_for_univariate_filter=1;
            T  = blkdiag(T,zeros(vobs));
            np    = size(T,1);
            Q   = blkdiag(Q,H);
            R  = blkdiag(R,eye(vobs));
            H   = zeros(vobs,vobs);
            Z   = [Z, eye(vobs)];
        end
    end
    [Pstar,Pinf] = compute_Pinf_Pstar(mf,T,R,Q,options_.qz_criterium);
elseif options_.lik_init == 4           % Start from the solution of the Riccati equation.
    [err, Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(mf,np,vobs)),H);
    mexErrCheck('kalman_steady_state',err);
    Pinf  = [];
    if kalman_algo~=2
        kalman_algo = 1;
    end
elseif options_.lik_init == 5            % Old diffuse Kalman filter only for the non stationary variables
    [eigenvect, eigenv] = eig(T);
    eigenv = diag(eigenv);
    nstable = length(find(abs(abs(eigenv)-1) > 1e-7));
    unstable = find(abs(abs(eigenv)-1) < 1e-7);
    V = eigenvect(:,unstable);
    indx_unstable = find(sum(abs(V),2)>1e-5);
    stable = find(sum(abs(V),2)<1e-5);
    nunit = length(eigenv) - nstable;
    Pstar = options_.Harvey_scale_factor*eye(np);
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    R_tmp = R(stable, :);
    T_tmp = T(stable,stable);
    Pstar_tmp=lyapunov_solver(T_tmp,R_tmp,Q,DynareOptions);
    Pstar(stable, stable) = Pstar_tmp;
    Pinf  = [];
end
kalman_tol = options_.kalman_tol;
diffuse_kalman_tol = options_.diffuse_kalman_tol;
riccati_tol = options_.riccati_tol;
data1 = Y-trend;
% -----------------------------------------------------------------------------
%  4. Kalman smoother
% -----------------------------------------------------------------------------

if ~missing_value
    for i=1:smpl
        data_index{i}=(1:vobs)';
    end
end

ST = T;
R1 = R;

if kalman_algo == 1 || kalman_algo == 3
    [alphahat,epsilonhat,etahat,ahat,P,aK,PK,decomp,state_uncertainty] = missing_DiffuseKalmanSmootherH1_Z(ST, ...
                                                      Z,R1,Q,H,Pinf,Pstar, ...
                                                      data1,vobs,np,smpl,data_index, ...
                                                      options_.nk,kalman_tol,diffuse_kalman_tol,options_.filter_decomposition,options_.smoothed_state_uncertainty);
    if isinf(alphahat)
        if kalman_algo == 1
            kalman_algo = 2;
        elseif kalman_algo == 3
            kalman_algo = 4;
        else
            error('This case shouldn''t happen')
        end
    end
end

if kalman_algo == 2 || kalman_algo == 4
    if ~all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
        if ~expanded_state_vector_for_univariate_filter
            %Augment state vector (follows Section 6.4.3 of DK (2012))
            expanded_state_vector_for_univariate_filter=1;
            Z   = [Z, eye(vobs)];
            ST  = blkdiag(ST,zeros(vobs));
            np  = size(ST,1);
            Q   = blkdiag(Q,H);
            R1  = blkdiag(R,eye(vobs));
            if kalman_algo == 4
                %recompute Schur state space transformation with
                %expanded state space
                [Pstar,Pinf] = compute_Pinf_Pstar(mf,ST,R1,Q,options_.qz_criterium);
            else
                Pstar = blkdiag(Pstar,H);
                if ~isempty(Pinf)
                    Pinf  = blkdiag(Pinf,zeros(vobs));
                end
            end
            %now reset H to 0
            H   = zeros(vobs,vobs);
        else
            %do nothing, state vector was already expanded
        end
    end

    [alphahat,epsilonhat,etahat,ahat,P,aK,PK,decomp,state_uncertainty] = missing_DiffuseKalmanSmootherH3_Z(ST, ...
                                                      Z,R1,Q,diag(H), ...
                                                      Pinf,Pstar,data1,vobs,np,smpl,data_index, ...
                                                      options_.nk,kalman_tol,diffuse_kalman_tol, ...
                                                      options_.filter_decomposition,options_.smoothed_state_uncertainty);
end


if expanded_state_vector_for_univariate_filter && (kalman_algo == 2 || kalman_algo == 4)
    % extracting measurement errors
    % removing observed variables from the state vector
    k = (1:np-vobs);
    alphahat = alphahat(k,:);
    ahat = ahat(k,:);
    aK = aK(:,k,:,:);
    epsilonhat=etahat(end-vobs+1:end,:);
    etahat=etahat(1:end-vobs,:);
    if ~isempty(PK)
        PK = PK(:,k,k,:);
    end
    if ~isempty(decomp)
        decomp = decomp(:,k,:,:);
    end
    if ~isempty(P)
        P = P(k,k,:);
    end
    if ~isempty(state_uncertainty)
        state_uncertainty = state_uncertainty(k,k,:);
    end
end
