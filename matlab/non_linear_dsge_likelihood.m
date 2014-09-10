function [fval,ys,trend_coeff,exit_flag,info,Model,DynareOptions,BayesInfo,DynareResults] = non_linear_dsge_likelihood(xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)
% Evaluates the posterior kernel of a dsge model using a non linear filter.

%@info:
%! @deftypefn {Function File} {[@var{fval},@var{exit_flag},@var{ys},@var{trend_coeff},@var{info},@var{Model},@var{DynareOptions},@var{BayesInfo},@var{DynareResults}] =} non_linear_dsge_likelihood (@var{xparam1},@var{DynareDataset},@var{DynareOptions},@var{Model},@var{EstimatedParameters},@var{BayesInfo},@var{DynareResults})
%! @anchor{dsge_likelihood}
%! @sp 1
%! Evaluates the posterior kernel of a dsge model using a non linear filter.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item xparam1
%! Vector of doubles, current values for the estimated parameters.
%! @item DynareDataset
%! Matlab's structure describing the dataset (initialized by dynare, see @ref{dataset_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item Model
%! Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%! @item EstimatedParamemeters
%! Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item fval
%! Double scalar, value of (minus) the likelihood.
%! @item exit_flag
%! Integer scalar, equal to zero if the routine return with a penalty (one otherwise).
%! @item ys
%! Vector of doubles, steady state level for the endogenous variables.
%! @item trend_coeffs
%! Matrix of doubles, coefficients of the deterministic trend in the measurement equation.
%! @item info
%! Integer scalar, error code.
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @item info==30
%! Ergodic variance can't be computed.
%! @item info==41
%! At least one parameter is violating a lower bound condition.
%! @item info==42
%! At least one parameter is violating an upper bound condition.
%! @item info==43
%! The covariance matrix of the structural innovations is not positive definite.
%! @item info==44
%! The covariance matrix of the measurement errors is not positive definite.
%! @item info==45
%! Likelihood is not a number (NaN).
%! @item info==45
%! Likelihood is a complex valued number.
%! @end table
%! @item Model
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dynare_estimation_1}, @ref{mode_check}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{dynare_resolve}, @ref{lyapunov_symm}, @ref{priordens}
%! @end deftypefn
%@eod:

% Copyright (C) 2010-2013 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr
%           frederic DOT karame AT univ DASH lemans DOT fr

global objective_function_penalty_base
% Declaration of the penalty as a persistent variable.
persistent init_flag
persistent restrict_variables_idx observed_variables_idx state_variables_idx mf0 mf1
persistent sample_size number_of_state_variables number_of_observed_variables number_of_structural_innovations

% Initialization of the returned arguments.
fval            = [];
ys              = [];
trend_coeff     = [];
exit_flag       = 1;

% Issue an error if loglinear option is used.
if DynareOptions.loglinear
    error('non_linear_dsge_likelihood: It is not possible to use a non linear filter with the option loglinear!')
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if (DynareOptions.mode_compute~=1) && any(xparam1<BayesInfo.lb)
    k = find(xparam1(:) < BayesInfo.lb);
    fval = objective_function_penalty_base+sum((BayesInfo.lb(k)-xparam1(k)).^2);
    exit_flag = 0;
    info = 41;
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if (DynareOptions.mode_compute~=1) && any(xparam1>BayesInfo.ub)
    k = find(xparam1(:)>BayesInfo.ub);
    fval = objective_function_penalty_base+sum((xparam1(k)-BayesInfo.ub(k)).^2);
    exit_flag = 0;
    info = 42;
    return
end

Model = set_all_parameters(xparam1,EstimatedParameters,Model);

Q = Model.Sigma_e;
H = Model.H;

if ~issquare(Q) || EstimatedParameters.ncx || isfield(EstimatedParameters,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Q);
    if ~Q_is_positive_definite
        fval = objective_function_penalty_base+penalty;
        exit_flag = 0;
        info = 43;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(EstimatedParameters.calibrated_covariances.position).^2);
            fval = objective_function_penalty_base+penalty;
            exit_flag = 0;
            info = 71;
            return
        end
    end

end

if ~issquare(H) || EstimatedParameters.ncn || isfield(EstimatedParameters,'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H);
    if ~H_is_positive_definite
        fval = objective_function_penalty_base+penalty;
        exit_flag = 0;
        info = 44;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances_ME')
        correct_flag=check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(EstimatedParameters.calibrated_covariances_ME.position).^2);
            fval = objective_function_penalty_base+penalty;
            exit_flag = 0;
            info = 72;
            return
        end
    end

end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic sdteadystate and extract the matrices of the state equation (T and R).
[T,R,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

if info(1) == 1 || info(1) == 2 || info(1) == 5 || info(1) == 25 || info(1) == 10
    fval = objective_function_penalty_base+1;
    exit_flag = 0;
    return
elseif info(1) == 3 || info(1) == 4 || info(1)==6 ||info(1) == 19 || info(1) == 20 || info(1) == 21
    fval = objective_function_penalty_base+info(2);
    exit_flag = 0;
    return
end

% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Define the deterministic linear trend of the measurement equation.
if DynareOptions.noconstant
    constant = zeros(DynareDataset.vobs,1);
else
    constant = SteadyState(BayesInfo.mfys);
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    trend_coeff = zeros(DynareDataset.vobs,1);
    t = DynareOptions.trend_coeffs;
    for i=1:length(t)
        if ~isempty(t{i})
            trend_coeff(i) = evalin('base',t{i});
        end
    end
    trend = repmat(constant,1,DynareDataset.nobs)+trend_coeff*[1:DynareDataset.nobs];
else
    trend = repmat(constant,1,DynareDataset.nobs);
end

% Get needed informations for kalman filter routines.
start = DynareOptions.presample+1;
np    = size(T,1);
mf    = BayesInfo.mf;
Y     = transpose(DynareDataset.data)-trend;

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------

% Get decision rules and transition equations.
dr = DynareResults.dr;

% Set persistent variables (first call).
if isempty(init_flag)
    mf0 = BayesInfo.mf0;
    mf1 = BayesInfo.mf1;
    restrict_variables_idx  = BayesInfo.restrict_var_list;
    observed_variables_idx  = restrict_variables_idx(mf1);
    state_variables_idx     = restrict_variables_idx(mf0);
    sample_size = size(Y,2);
    number_of_state_variables = length(mf0);
    number_of_observed_variables = length(mf1);
    number_of_structural_innovations = length(Q);
    init_flag = 1;
end

ReducedForm.ghx  = dr.ghx(restrict_variables_idx,:);
ReducedForm.ghu  = dr.ghu(restrict_variables_idx,:);
ReducedForm.ghxx = dr.ghxx(restrict_variables_idx,:);
ReducedForm.ghuu = dr.ghuu(restrict_variables_idx,:);
ReducedForm.ghxu = dr.ghxu(restrict_variables_idx,:);
ReducedForm.steadystate = dr.ys(dr.order_var(restrict_variables_idx));
ReducedForm.constant = ReducedForm.steadystate + .5*dr.ghs2(restrict_variables_idx);
ReducedForm.state_variables_steady_state = dr.ys(dr.order_var(state_variables_idx));
ReducedForm.Q = Q;
ReducedForm.H = H;
ReducedForm.mf0 = mf0;
ReducedForm.mf1 = mf1;

% Set initial condition.
switch DynareOptions.particle.initialization
  case 1% Initial state vector covariance is the ergodic variance associated to the first order Taylor-approximation of the model.
    StateVectorMean = ReducedForm.constant(mf0);
    StateVectorVariance = lyapunov_symm(ReducedForm.ghx(mf0,:),ReducedForm.ghu(mf0,:)*ReducedForm.Q*ReducedForm.ghu(mf0,:)',1e-12,1e-12);
  case 2% Initial state vector covariance is a monte-carlo based estimate of the ergodic variance (consistent with a k-order Taylor-approximation of the model).
    StateVectorMean = ReducedForm.constant(mf0);
    old_DynareOptionsperiods = DynareOptions.periods;
    DynareOptions.periods = 5000;
    y_ = simult(oo_.steady_state, dr,Model,DynareOptions,DynareResults);
    y_ = y_(state_variables_idx,2001:5000);
    StateVectorVariance = cov(y_');
    DynareOptions.periods = old_DynareOptionsperiods;
    clear('old_DynareOptionsperiods','y_');
  case 3% Initial state vector covariance is a diagonal matrix.
    StateVectorMean = ReducedForm.constant(mf0);
    StateVectorVariance = DynareOptions.particle.initial_state_prior_std*eye(number_of_state_variables);
  otherwise
    error('Unknown initialization option!')
end
ReducedForm.StateVectorMean = StateVectorMean;
ReducedForm.StateVectorVariance = StateVectorVariance;

%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------
DynareOptions.warning_for_steadystate = 0;
[s1,s2] = get_dynare_random_generator_state();
LIK = feval(DynareOptions.particle.algorithm,ReducedForm,Y,start,DynareOptions);
set_dynare_random_generator_state(s1,s2);
if imag(LIK)
    likelihood = objective_function_penalty_base;
    exit_flag  = 0;
elseif isnan(LIK)
    likelihood = objective_function_penalty_base;
    exit_flag  = 0;
else
    likelihood = LIK;
end
DynareOptions.warning_for_steadystate = 1;
% ------------------------------------------------------------------------------
% Adds prior if necessary
% ------------------------------------------------------------------------------
lnprior = priordens(xparam1(:),BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
fval    = (likelihood-lnprior);