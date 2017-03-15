function [oo_,M_,options_,bayestopt_,Smoothed_variables_declaration_order_deviation_form]=evaluate_smoother(parameters,var_list,M_,oo_,options_,bayestopt_,estim_params_)
% Evaluate the smoother at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean','mle_mode') or a vector of values for
%                  the (estimated) parameters of the model.
%    o var_list    subset of endogenous variables
%    o M_          [structure]  Definition of the model
%    o oo_         [structure]  Storage of results
%    o options_    [structure]  Options
%    o bayestopt_  [structure]  describing the priors
%    o estim_params_ [structure] characterizing parameters to be estimated
%
% OUTPUTS
%    o oo       [structure]  results:
%                              - SmoothedVariables
%                              - SmoothedShocks
%                              - FilteredVariablesShockDecomposition
%                              - UpdatedVariables
%                              - FilteredVariables
%                              - SmoothedMeasurementErrors
%                              - FilteredVariablesKStepAhead
%                              - FilteredVariablesKStepAheadVariances
%    o M_          [structure]  Definition of the model
%    o options_    [structure]  Options; returns options_.first_obs
%    o bayestopt_  [structure]  describing the priors; returns fields like bayestopt_.smoother_var_list from the smoother 
%    o Smoothed_variables_declaration_order_deviation_form
%                           Smoothed variables from the Kalman smoother in
%                           order of declaration of variables (M_.endo_names)
%                           in deviations from their respective mean, i.e.
%                           without trend and constant part (used for shock_decomposition)
%
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function use persistent variables for the dataset and the description of the missing observations. Consequently, if this function
%     is called more than once (by changing the value of parameters) the sample *must not* change.

% Copyright (C) 2010-2016 Dynare Team
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

persistent dataset_ dataset_info

%store qz_criterium
qz_criterium_old=options_.qz_criterium;

if ischar(parameters) && strcmp(parameters,'calibration')
    options_.smoother=1;
end

if isempty(dataset_) || isempty(bayestopt_)
    [dataset_,dataset_info,xparam1, hh, M_, options_, oo_, estim_params_,bayestopt_] = dynare_estimation_init(var_list, M_.fname, [], M_, options_, oo_, estim_params_, bayestopt_);
end

if nargin==0
    parameters = 'posterior_mode';
end

if ischar(parameters)
    switch parameters
      case 'posterior_mode'
        parameters = get_posterior_parameters('mode');
      case 'posterior_mean'
        parameters = get_posterior_parameters('mean');
      case 'posterior_median'
        parameters = get_posterior_parameters('median');
      case 'mle_mode'
        parameters = get_posterior_parameters('mode','mle_');
      case 'prior_mode'
        parameters = bayestopt_.p5(:);
      case 'prior_mean'
        parameters = bayestopt_.p1;
      case 'calibration'
        if isempty(oo_.dr)
            error('You must run ''stoch_simul'' first.');
        end
        parameters = [];
      otherwise
        disp('evaluate_smoother:: If the input argument is a string, then it has to be equal to:')
        disp('                     ''posterior_mode'', ')
        disp('                     ''posterior_mean'', ')
        disp('                     ''posterior_median'', ')
        disp('                     ''prior_mode'' or')
        disp('                     ''prior_mean''.')
        disp('                     ''calibration''.')
        error
    end
end

[atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,M_,oo_,options_,bayestopt_] = ...
    DsgeSmoother(parameters,dataset_.nobs,transpose(dataset_.data),dataset_info.missing.aindex,dataset_info.missing.state,M_,oo_,options_,bayestopt_,estim_params_);
[oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);

if nargout==4
   Smoothed_variables_declaration_order_deviation_form=atT(oo_.dr.inv_order_var(bayestopt_.smoother_var_list),:);
end

%reset qz_criterium
options_.qz_criterium=qz_criterium_old;
