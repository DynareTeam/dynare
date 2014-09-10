function [xparams,lpd,hessian] = ...
    maximize_prior_density(iparams, prior_shape, prior_hyperparameter_1, prior_hyperparameter_2, prior_inf_bound, prior_sup_bound,DynareOptions,DynareModel,EstimatedParams,DynareResults)
% Maximizes the logged prior density using Chris Sims' optimization routine.
% 
% INPUTS 
%   iparams                        [double]   vector of initial parameters.
%   prior_shape                    [integer]  vector specifying prior densities shapes.
%   prior_hyperparameter_1         [double]   vector, first hyperparameter.
%   prior_hyperparameter_2         [double]   vector, second hyperparameter.
%   prior_inf_bound                [double]   vector, prior's lower bound.
%   prior_sup_bound                [double]   vector, prior's upper bound. 
%
% OUTPUTS 
%   xparams       [double]  vector, prior mode.
%   lpd           [double]  scalar, value of the logged prior density at the mode.
%   hessian       [double]  matrix, Hessian matrix at the prior mode.

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

number_of_estimated_parameters = length(iparams);
H0 = 1e-4*eye(number_of_estimated_parameters);
crit = 1e-7;
nit = 1000;
gradient_method = 2;
gradient_epsilon = 1e-6;

[lpd,xparams,grad,hessian,itct,fcount,retcodehat] = ...
    csminwel1('minus_logged_prior_density',iparams,H0,[],crit,nit,gradient_method, gradient_epsilon, ...
              prior_shape, prior_hyperparameter_1, prior_hyperparameter_2, prior_inf_bound, prior_sup_bound,DynareOptions,DynareModel,EstimatedParams,DynareResults);

lpd = -lpd;
