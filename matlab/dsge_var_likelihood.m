function [fval,grad,hess,exit_flag,info,PHI,SIGMAu,iXX,prior] = dsge_var_likelihood(xparam1,DynareDataset,DynareInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)
% Evaluates the posterior kernel of the bvar-dsge model.
%
% INPUTS
%   o xparam1       [double]     Vector of model's parameters.
%   o gend          [integer]    Number of observations (without conditionning observations for the lags).
%
% OUTPUTS
%   o fval          [double]     Value of the posterior kernel at xparam1.
%   o cost_flag     [integer]    Zero if the function returns a penalty, one otherwise.
%   o info          [integer]    Vector of informations about the penalty.
%   o PHI           [double]     Stacked BVAR-DSGE autoregressive matrices (at the mode associated to xparam1).
%   o SIGMAu        [double]     Covariance matrix of the BVAR-DSGE (at the mode associated to xparam1).
%   o iXX           [double]     inv(X'X).
%   o prior         [double]     a matlab structure describing the dsge-var prior.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2012 Dynare Team
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

global objective_function_penalty_base
persistent dsge_prior_weight_idx

grad=[];
hess=[];
exit_flag = [];
info = [];
PHI = [];
SIGMAu = [];
iXX = [];
prior = [];

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

% Initialize some of the output arguments.
fval = [];
exit_flag = 1;

% Return, with endogenous penalty, if some dsge-parameters are smaller than the lower bound of the prior domain.
if DynareOptions.mode_compute ~= 1 && any(xparam1 < BayesInfo.lb)
    k = find(xparam1 < BayesInfo.lb);
    fval = objective_function_penalty_base+sum((BayesInfo.lb(k)-xparam1(k)).^2);
    exit_flag = 0;
    info = 41;
    return;
end

% Return, with endogenous penalty, if some dsge-parameters are greater than the upper bound of the prior domain.
if DynareOptions.mode_compute ~= 1 && any(xparam1 > BayesInfo.ub)
    k = find(xparam1 > BayesInfo.ub);
    fval = objective_function_penalty_base+sum((xparam1(k)-BayesInfo.ub(k)).^2);
    exit_flag = 0;
    info = 42;
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
    fval = objective_function_penalty_base+abs(NumberOfObservations*dsge_prior_weight-(NumberOfParameters+NumberOfObservedVariables));
    exit_flag = 0;
    info = 51;
    info(2)=dsge_prior_weight;
    info(3)=(NumberOfParameters+NumberOfObservedVariables)/DynareDataset.nobs;
    return
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Solve the Dsge model and get the matrices of the reduced form solution. T and R are the matrices of the
% state equation
[T,R,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1) == 1 || info(1) == 2 || info(1) == 5 || info(1) == 7 || info(1) ...
            == 8 || info(1) == 22 || info(1) == 24 || info(1) == 25 ...
            || info(1) == 10
    fval = objective_function_penalty_base+1;
    info = info(1);
    exit_flag = 0;
    return
elseif info(1) == 3 || info(1) == 4 || info(1) == 19 || info(1) == 20 || info(1) == 21
    fval = objective_function_penalty_base+info(2);
    info = info(1);
    exit_flag = 0;
    return
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
tmp0 = lyapunov_symm(T,R*Q*R',DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
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

if ~isinf(dsge_prior_weight)% Evaluation of the likelihood of the dsge-var model when the dsge prior weight is finite.
    tmp0 = dsge_prior_weight*NumberOfObservations*TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1) + mYY ;
    tmp1 = dsge_prior_weight*NumberOfObservations*GYX + mYX;
    tmp2 = inv(dsge_prior_weight*NumberOfObservations*GXX+mXX);
    SIGMAu = tmp0 - tmp1*tmp2*tmp1'; clear('tmp0');
    [SIGMAu_is_positive_definite, penalty] = ispd(SIGMAu);
    if ~SIGMAu_is_positive_definite
        fval = objective_function_penalty_base + penalty;
        info = 52;
        exit_flag = 0;
        return;
    end
    SIGMAu = SIGMAu / (NumberOfObservations*(1+dsge_prior_weight));
    PHI = tmp2*tmp1'; clear('tmp1');
    prodlng1 = sum(gammaln(.5*((1+dsge_prior_weight)*NumberOfObservations- ...
                               NumberOfObservedVariables*NumberOfLags ...
                               +1-(1:NumberOfObservedVariables)')));
    prodlng2 = sum(gammaln(.5*(dsge_prior_weight*NumberOfObservations- ...
                               NumberOfObservedVariables*NumberOfLags ...
                               +1-(1:NumberOfObservedVariables)')));
    lik = .5*NumberOfObservedVariables*log(det(dsge_prior_weight*NumberOfObservations*GXX+mXX)) ...
          + .5*((dsge_prior_weight+1)*NumberOfObservations-NumberOfParameters)*log(det((dsge_prior_weight+1)*NumberOfObservations*SIGMAu)) ...
          - .5*NumberOfObservedVariables*log(det(dsge_prior_weight*NumberOfObservations*GXX)) ...
          - .5*(dsge_prior_weight*NumberOfObservations-NumberOfParameters)*log(det(dsge_prior_weight*NumberOfObservations*(GYY-GYX*inv(GXX)*GYX'))) ...
          + .5*NumberOfObservedVariables*NumberOfObservations*log(2*pi)  ...
          - .5*log(2)*NumberOfObservedVariables*((dsge_prior_weight+1)*NumberOfObservations-NumberOfParameters) ...
          + .5*log(2)*NumberOfObservedVariables*(dsge_prior_weight*NumberOfObservations-NumberOfParameters) ...
          - prodlng1 + prodlng2;
else% Evaluation of the likelihood of the dsge-var model when the dsge prior weight is infinite.
    iGXX = inv(GXX);
    SIGMAu = GYY - GYX*iGXX*transpose(GYX);
    PHI = iGXX*transpose(GYX);
    lik = NumberOfObservations * ( log(det(SIGMAu)) + NumberOfObservedVariables*log(2*pi) +  ...
                   trace(inv(SIGMAu)*(mYY - transpose(mYX*PHI) - mYX*PHI + transpose(PHI)*mXX*PHI)/NumberOfObservations));
    lik = .5*lik;% Minus likelihood
end

% Add the (logged) prior density for the dsge-parameters.
lnprior = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
fval = (lik-lnprior);

if (nargout == 8)
    if isinf(dsge_prior_weight)
        iXX = iGXX;
    else
        iXX = tmp2;
    end
end

if (nargout==9)
    if isinf(dsge_prior_weight)
        iXX = iGXX;
    else
        iXX = tmp2;
    end
    iGXX = inv(GXX);
    prior.SIGMAstar = GYY - GYX*iGXX*GYX';
    prior.PHIstar = iGXX*transpose(GYX);
    prior.ArtificialSampleSize = fix(dsge_prior_weight*NumberOfObservations);
    prior.DF = prior.ArtificialSampleSize - NumberOfParameters - NumberOfObservedVariables;
    prior.iGXX = iGXX;
end