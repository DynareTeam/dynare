function [fval,fake_1, fake_2, exit_flag ] = minus_logged_prior_density(xparams,pshape,p6,p7,p3,p4,DynareOptions,DynareModel,EstimatedParams,DynareResults)
% Evaluates minus the logged prior density.
% 
% INPUTS 
%   xparams    [double]   vector of parameters.
%   pshape     [integer]  vector specifying prior densities shapes.
%   p6         [double]   vector, first hyperparameter.
%   p7         [double]   vector, second hyperparameter.
%   p3         [double]   vector, prior's lower bound.
%   p4         [double]   vector, prior's upper bound. 
%
% OUTPUTS 
%   f          [double]  value of minus the logged prior density.

% Copyright (C) 2009-2012 Dynare Team
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

fake_1 = 1;
fake_2 = 1;

exit_flag = 1;
info = 0;

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if ~isequal(DynareOptions.mode_compute,1) && any(xparams<p3)
    k = find(xparams<p3);
    fval = objective_function_penalty_base+sum((p3(k)-xparams(k)).^2);
    exit_flag = 0;
    info = 41;
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if ~isequal(DynareOptions.mode_compute,1) && any(xparams>p4)
    k = find(xparams>p4);
    fval = objective_function_penalty_base+sum((xparams(k)-p4(k)).^2);
    exit_flag = 0;
    info = 42;
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
DynareModel = set_all_parameters(xparams,EstimatedParams,DynareModel);

Q = DynareModel.Sigma_e;
H = DynareModel.H;

% Test if Q is positive definite.
if ~issquare(Q) || EstimatedParams.ncx || isfield(EstimatedParams,'calibrated_covariances')
    % Try to compute the cholesky decomposition of Q (possible iff Q is positive definite)
    [Q_is_positive_definite, penalty] = ispd(Q);
    if ~Q_is_positive_definite
        % The variance-covariance matrix of the structural innovations is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
        fval = objective_function_penalty_base+penalty;
        exit_flag = 0;
        info = 43;
        return
    end
    if isfield(EstimatedParams,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(EstimatedParams.calibrated_covariances.position).^2);
            fval = objective_function_penalty_base+penalty;
            exit_flag = 0;
            info = 71;
            return
        end
    end

end

% Test if H is positive definite.
if ~issquare(H) || EstimatedParams.ncn || isfield(EstimatedParams,'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H);
    if ~H_is_positive_definite
        % The variance-covariance matrix of the measurement errors is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
        fval = objective_function_penalty_base+penalty;
        exit_flag = 0;
        info = 44;
        return
    end
    if isfield(EstimatedParams,'calibrated_covariances_ME')
        correct_flag=check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(EstimatedParams.calibrated_covariances_ME.position).^2);
            fval = objective_function_penalty_base+penalty;
            exit_flag = 0;
            info = 72;
            return
        end
    end
end


%-----------------------------
% 2. Check BK and steady state
%-----------------------------

M_ = set_all_parameters(xparams,EstimatedParams,DynareModel);
[dr,info,DynareModel,DynareOptions,DynareResults] = resol(0,DynareModel,DynareOptions,DynareResults);

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1) == 1 || info(1) == 2 || info(1) == 5 || info(1) == 7 || info(1) ...
            == 8 || info(1) == 22 || info(1) == 24 || info(1) == 19
    fval = objective_function_penalty_base+1;
    info = info(1);
    exit_flag = 0;
    return
elseif info(1) == 3 || info(1) == 4 || info(1)==6 || info(1) == 20 || info(1) == 21  || info(1) == 23
    fval = objective_function_penalty_base+info(2);
    info = info(1);
    exit_flag = 0;
    return
end


fval = - priordens(xparams,pshape,p6,p7,p3,p4);