function [fval,info,exit_flag,grad,hess,SteadyState,trend_coeff,PHI_tilde,SIGMA_u_tilde,iXX,prior] = dsge_var_likelihood(xparam1,DynareDataset,DynareInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults)
% Evaluates the posterior kernel of the bvar-dsge model.
%
% INPUTS
%   o xparam1       [double]     Vector of model's parameters.
%   o gend          [integer]    Number of observations (without conditionning observations for the lags).
%
% OUTPUTS
%   o fval          [double]     Value of the posterior kernel at xparam1.
%   o info          [integer]    Vector of informations about the penalty.
%   o exit_flag     [integer]    Zero if the function returns a penalty, one otherwise.
%   o grad          [double]     place holder for gradient of the likelihood 
%                                currently not supported by dsge_var
%   o hess          [double]     place holder for hessian matrix of the likelihood 
%                                currently not supported by dsge_var
%   o SteadyState   [double]     Steady state vector possibly recomputed
%                                by call to dynare_resolve()
%   o trend_coeff   [double]     place holder for trend coefficients,
%                                currently not supported by dsge_var
%   o PHI_tilde     [double]     Stacked BVAR-DSGE autoregressive matrices (at the mode associated to xparam1); 
%                                formula (28), DS (2004)
%   o SIGMA_u_tilde [double]     Covariance matrix of the BVAR-DSGE (at the mode associated to xparam1),
%                                formula (29), DS (2004)
%   o iXX           [double]     inv(lambda*T*Gamma_XX^*+ X'*X)
%   o prior         [double]     a matlab structure describing the dsge-var prior
%                                   - SIGMA_u_star: prior covariance matrix, formula (23), DS (2004)
%                                   - PHI_star: prior autoregressive matrices, formula (22), DS (2004)
%                                   - ArtificialSampleSize: number of artificial observations from the prior (T^* in DS (2004))
%                                   - DF = prior.ArtificialSampleSize - NumberOfParameters - NumberOfObservedVariables;
%                                   - iGXX_star: theoretical covariance of regressors ({\Gamma_{XX}^*}^{-1} in DS (2004))
%
% ALGORITHMS
%   Follows the computations outlined in Del Negro/Schorfheide (2004):
%   Priors from general equilibrium models for VARs, International Economic
%   Review, 45(2), pp. 643-673 
% 
% SPECIAL REQUIREMENTS
%   None.

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

persistent dsge_prior_weight_idx

% Initialize some of the output arguments.
fval = [];
exit_flag = 1;
grad=[];
hess=[];
info = zeros(4,1);
PHI_tilde = [];
SIGMA_u_tilde = [];
iXX = [];
prior = [];
trend_coeff=[];

% Ensure that xparam1 is a column vector.
xparam1 = xparam1(:);

% Initialization of of the index for parameter dsge_prior_weight in Model.params.
if isempty(dsge_prior_weight_idx)
    dsge_prior_weight_idx = strmatch('dsge_prior_weight',Model.param_names);
end

% Get the number of estimated (dsge) parameters.
nx = EstimatedParameters.nvx + EstimatedParameters.np;

% Get the number of observed variables in the VAR model.
NumberOfObservedVariables = DynareDataset.vobs;

% Get the number of observations.
NumberOfObservations = DynareDataset.nobs;


% Get the number of lags in the VAR model.
NumberOfLags = DynareOptions.dsge_varlag;

% Get the number of parameters in the VAR model.
NumberOfParameters = NumberOfObservedVariables*NumberOfLags ;
if ~DynareOptions.noconstant
    NumberOfParameters = NumberOfParameters + 1;
end

% Get empirical second order moments for the observed variables.
mYY = evalin('base', 'mYY');
mYX = evalin('base', 'mYX');
mXY = evalin('base', 'mXY');
mXX = evalin('base', 'mXX');

% Return, with endogenous penalty, if some dsge-parameters are smaller than the lower bound of the prior domain.
if DynareOptions.mode_compute ~= 1 && any(xparam1 < BoundsInfo.lb)
    k = find(xparam1 < BoundsInfo.lb);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4)= sum((BoundsInfo.lb(k)-xparam1(k)).^2);
    return;
end

% Return, with endogenous penalty, if some dsge-parameters are greater than the upper bound of the prior domain.
if DynareOptions.mode_compute ~= 1 && any(xparam1 > BoundsInfo.ub)
    k = find(xparam1 > BoundsInfo.ub);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4) = sum((xparam1(k)-BoundsInfo.ub(k)).^2);
    return;
end

% Get the variance of each structural innovation.
Q = Model.Sigma_e;
for i=1:EstimatedParameters.nvx
    k = EstimatedParameters.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
end
offset = EstimatedParameters.nvx;

% Update Model.params and Model.Sigma_e.
Model.params(EstimatedParameters.param_vals(:,1)) = xparam1(offset+1:end);
Model.Sigma_e = Q;

% Get the weight of the dsge prior.
dsge_prior_weight = Model.params(dsge_prior_weight_idx);

% Is the dsge prior proper?
if dsge_prior_weight<(NumberOfParameters+NumberOfObservedVariables)/NumberOfObservations;
    fval = Inf;
    exit_flag = 0;
    info(1) = 51;
    info(2)=dsge_prior_weight;
    info(3)=(NumberOfParameters+NumberOfObservedVariables)/DynareDataset.nobs;
    info(4)=abs(NumberOfObservations*dsge_prior_weight-(NumberOfParameters+NumberOfObservedVariables));
    return
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Solve the Dsge model and get the matrices of the reduced form solution. T and R are the matrices of the
% state equation
[T,R,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
            info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
            info(1) == 81 || info(1) == 84 ||  info(1) == 85
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end
end

% Define the mean/steady state vector.
if ~DynareOptions.noconstant
    if DynareOptions.loglinear
        constant = transpose(log(SteadyState(BayesInfo.mfys)));
    else
        constant = transpose(SteadyState(BayesInfo.mfys));
    end
else
    constant = zeros(1,NumberOfObservedVariables);
end


%------------------------------------------------------------------------------
% 3. theoretical moments (second order)
%------------------------------------------------------------------------------

% Compute the theoretical second order moments
tmp0 = lyapunov_symm(T,R*Q*R',DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold, [], DynareOptions.debug);
mf  = BayesInfo.mf1;

% Get the non centered second order moments
TheoreticalAutoCovarianceOfTheObservedVariables = zeros(NumberOfObservedVariables,NumberOfObservedVariables,NumberOfLags+1);
TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1) = tmp0(mf,mf)+constant'*constant;
for lag = 1:NumberOfLags
    tmp0 = T*tmp0;
    TheoreticalAutoCovarianceOfTheObservedVariables(:,:,lag+1) = tmp0(mf,mf) + constant'*constant;
end

% Build the theoretical "covariance" between Y and X
GYX = zeros(NumberOfObservedVariables,NumberOfParameters);
for i=1:NumberOfLags
    GYX(:,(i-1)*NumberOfObservedVariables+1:i*NumberOfObservedVariables) = TheoreticalAutoCovarianceOfTheObservedVariables(:,:,i+1);
end
if ~DynareOptions.noconstant
    GYX(:,end) = constant';
end

% Build the theoretical "covariance" between X and X
GXX = kron(eye(NumberOfLags), TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1));
for i = 1:NumberOfLags-1
    tmp1 = diag(ones(NumberOfLags-i,1),i);
    tmp2 = diag(ones(NumberOfLags-i,1),-i);
    GXX = GXX + kron(tmp1,TheoreticalAutoCovarianceOfTheObservedVariables(:,:,i+1));
    GXX = GXX + kron(tmp2,TheoreticalAutoCovarianceOfTheObservedVariables(:,:,i+1)');
end

if ~DynareOptions.noconstant
    % Add one row and one column to GXX
    GXX = [GXX , kron(ones(NumberOfLags,1),constant') ; [  kron(ones(1,NumberOfLags),constant) , 1] ];
end

GYY = TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1);

assignin('base','GYY',GYY);
assignin('base','GXX',GXX);
assignin('base','GYX',GYX);

iGXX = inv(GXX);
PHI_star = iGXX*transpose(GYX); %formula (22), DS (2004)
SIGMA_u_star=GYY - GYX*PHI_star; %formula (23), DS (2004)
if ~isinf(dsge_prior_weight)% Evaluation of the likelihood of the dsge-var model when the dsge prior weight is finite.
    tmp0 = dsge_prior_weight*NumberOfObservations*TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1) + mYY ;  %first term of square bracket in formula (29), DS (2004)
    tmp1 = dsge_prior_weight*NumberOfObservations*GYX + mYX;        %first element of second term of square bracket in formula (29), DS (2004)
    tmp2 = inv(dsge_prior_weight*NumberOfObservations*GXX+mXX);     %middle element of second term of square bracket in formula (29), DS (2004)
    SIGMA_u_tilde = tmp0 - tmp1*tmp2*tmp1';                               %square bracket term in formula (29), DS (2004) 
    clear('tmp0');
    [SIGMAu_is_positive_definite, penalty] = ispd(SIGMA_u_tilde);
    if ~SIGMAu_is_positive_definite
        fval = Inf;
        info(1) = 52;
        info(4) = penalty;
        exit_flag = 0;
        return;
    end
    SIGMA_u_tilde = SIGMA_u_tilde / (NumberOfObservations*(1+dsge_prior_weight));   %prefactor of formula (29), DS (2004)
    PHI_tilde = tmp2*tmp1';                                                   %formula (28), DS (2004)
    clear('tmp1');
    prodlng1 = sum(gammaln(.5*((1+dsge_prior_weight)*NumberOfObservations- ...
                               NumberOfParameters ...
                               +1-(1:NumberOfObservedVariables)')));    %last term in numerator of third line of (A.2), DS (2004)
    prodlng2 = sum(gammaln(.5*(dsge_prior_weight*NumberOfObservations- ...
                               NumberOfParameters ...
                               +1-(1:NumberOfObservedVariables)')));    %last term in denominator of third line of (A.2), DS (2004)
    %Compute minus log likelihood according to (A.2), DS (2004)
    lik = .5*NumberOfObservedVariables*log(det(dsge_prior_weight*NumberOfObservations*GXX+mXX)) ... %first term in numerator of second line of (A.2), DS (2004)
          + .5*((dsge_prior_weight+1)*NumberOfObservations-NumberOfParameters)*log(det((dsge_prior_weight+1)*NumberOfObservations*SIGMA_u_tilde)) ... %second term in numerator of second line of (A.2), DS (2004)
          - .5*NumberOfObservedVariables*log(det(dsge_prior_weight*NumberOfObservations*GXX)) ... %first term in denominator of second line of (A.2), DS (2004)
          - .5*(dsge_prior_weight*NumberOfObservations-NumberOfParameters)*log(det(dsge_prior_weight*NumberOfObservations*SIGMA_u_star)) ... %second term in denominator of second line of (A.2), DS (2004)
          + .5*NumberOfObservedVariables*NumberOfObservations*log(2*pi)  ... %first term in numerator of third line of (A.2), DS (2004)
          - .5*log(2)*NumberOfObservedVariables*((dsge_prior_weight+1)*NumberOfObservations-NumberOfParameters) ... %second term in numerator of third line of (A.2), DS (2004)
          + .5*log(2)*NumberOfObservedVariables*(dsge_prior_weight*NumberOfObservations-NumberOfParameters) ... %first term in denominator of third line of (A.2), DS (2004)
          - prodlng1 + prodlng2;
else% Evaluation of the likelihood of the dsge-var model when the dsge prior weight is infinite.
    PHI_star = iGXX*transpose(GYX);
    %Compute minus log likelihood according to (33), DS (2004) (where the last term in the trace operator has been multiplied out)
    lik = NumberOfObservations * ( log(det(SIGMA_u_star)) + NumberOfObservedVariables*log(2*pi) +  ...
                   trace(inv(SIGMA_u_star)*(mYY - transpose(mYX*PHI_star) - mYX*PHI_star + transpose(PHI_star)*mXX*PHI_star)/NumberOfObservations));
    lik = .5*lik;% Minus likelihood
    SIGMA_u_tilde=SIGMA_u_star;
    PHI_tilde=PHI_star;
end

if isnan(lik)
    fval = Inf;
    info(1) = 45;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if imag(lik)~=0
    fval = Inf;
    info(1) = 46;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

% Add the (logged) prior density for the dsge-parameters.
lnprior = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
fval = (lik-lnprior);

if isnan(fval)
    fval = Inf;
    info(1) = 47;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if imag(fval)~=0
    fval = Inf;
    info(1) = 48;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if isinf(fval)~=0
    fval = Inf;
    info(1) = 50;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if (nargout >= 10)
    if isinf(dsge_prior_weight)
        iXX = iGXX;
    else
        iXX = tmp2;
    end
end

if (nargout==11)
    prior.SIGMA_u_star = SIGMA_u_star; 
    prior.PHI_star = PHI_star; 
    prior.ArtificialSampleSize = fix(dsge_prior_weight*NumberOfObservations);
    prior.DF = prior.ArtificialSampleSize - NumberOfParameters - NumberOfObservedVariables;
    prior.iGXX_star = iGXX;
end