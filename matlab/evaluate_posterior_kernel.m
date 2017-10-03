function lpkern = evaluate_posterior_kernel(parameters,M_,estim_params_,oo_,options_,bayestopt_,llik)
% Evaluate the evaluate_posterior_kernel at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean') or a vector of values for
%                  the (estimated) parameters of the model.
%    o M_         [structure]     Definition of the model
%    o estim_params_ [structure]  characterizing parameters to be estimated
%    o oo_         [structure]    Storage of results
%    o options_    [structure]    Options
%    o bayestopt_  [structure]    describing the priors
%    o llik        [double]       value of the logged likelihood if it
%                                   should not be computed
%
% OUTPUTS
%    o lpkern      [double]  value of the logged posterior kernel.
%
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function cannot evaluate the prior density of a dsge-var model...
% [2] This function use persistent variables for the dataset and the description of the missing observations. Consequently, if this function
%     is called more than once (by changing the value of parameters) the sample *must not* change.

% Copyright (C) 2009-2017 Dynare Team
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

[ldens,parameters] = evaluate_prior(parameters,M_,estim_params_,oo_,options_,bayestopt_);
if nargin==6 %llik provided as an input
    llik = evaluate_likelihood(parameters,M_,estim_params_,oo_,options_,bayestopt_);
end
lpkern = ldens+llik;