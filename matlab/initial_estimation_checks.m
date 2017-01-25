function DynareResults = initial_estimation_checks(objective_function,xparam1,DynareDataset,DatasetInfo,Model,EstimatedParameters,DynareOptions,BayesInfo,BoundsInfo,DynareResults)
% function DynareResults = initial_estimation_checks(objective_function,xparam1,DynareDataset,DatasetInfo,Model,EstimatedParameters,DynareOptions,BayesInfo,BoundsInfo,DynareResults)
% Checks data (complex values, ML evaluation, initial values, BK conditions,..)
%
% INPUTS
%   objective_function  [function handle] of the objective function
%   xparam1:            [vector] of parameters to be estimated
%   DynareDataset:      [dseries] object storing the dataset
%   DataSetInfo:        [structure] storing informations about the sample.
%   Model:              [structure] decribing the model
%   EstimatedParameters [structure] characterizing parameters to be estimated
%   DynareOptions       [structure] describing the options
%   BayesInfo           [structure] describing the priors
%   BoundsInfo          [structure] containing prior bounds
%   DynareResults       [structure] storing the results
%
% OUTPUTS
%    DynareResults     structure of temporary results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2016 Dynare Team
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

%get maximum number of simultaneously observed variables for stochastic
%singularity check
maximum_number_non_missing_observations=max(sum(~isnan(DynareDataset.data),2));

if DynareOptions.order>1 && any(any(isnan(DynareDataset.data)))
    error('initial_estimation_checks:: particle filtering does not support missing observations')
end

if maximum_number_non_missing_observations>Model.exo_nbr+EstimatedParameters.nvn
    error(['initial_estimation_checks:: Estimation can''t take place because there are less declared shocks than observed variables!'])
end

if maximum_number_non_missing_observations>length(find(diag(Model.Sigma_e)))+EstimatedParameters.nvn
    error(['initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance!'])
end

if (any(BayesInfo.pshape  >0 ) && DynareOptions.mh_replic) && DynareOptions.mh_nblck<1
    error(['initial_estimation_checks:: Bayesian estimation cannot be conducted with mh_nblocks=0.'])    
end

old_steady_params=Model.params; %save initial parameters for check if steady state changes param values

% % check if steady state solves static model (except if diffuse_filter == 1)
[DynareResults.steady_state, new_steady_params] = evaluate_steady_state(DynareResults.steady_state,Model,DynareOptions,DynareResults,DynareOptions.diffuse_filter==0);

if isfield(EstimatedParameters,'param_vals') && ~isempty(EstimatedParameters.param_vals)
    %check whether steady state file changes estimated parameters
    Model_par_varied=Model; %store Model structure
    Model_par_varied.params(EstimatedParameters.param_vals(:,1))=Model_par_varied.params(EstimatedParameters.param_vals(:,1))*1.01; %vary parameters
    [junk, new_steady_params_2] = evaluate_steady_state(DynareResults.steady_state,Model_par_varied,DynareOptions,DynareResults,DynareOptions.diffuse_filter==0);

    changed_par_indices=find((old_steady_params(EstimatedParameters.param_vals(:,1))-new_steady_params(EstimatedParameters.param_vals(:,1))) ...
            | (Model_par_varied.params(EstimatedParameters.param_vals(:,1))-new_steady_params_2(EstimatedParameters.param_vals(:,1))));

    if ~isempty(changed_par_indices)
        fprintf('\nThe steady state file internally changed the values of the following estimated parameters:\n')
        disp(Model.param_names(EstimatedParameters.param_vals(changed_par_indices,1),:));
        fprintf('This will override the parameter values drawn from the proposal density and may lead to wrong results.\n')
        fprintf('Check whether this is really intended.\n')    
        warning('The steady state file internally changes the values of the estimated parameters.')
    end
end

if any(BayesInfo.pshape) % if Bayesian estimation
    nvx=EstimatedParameters.nvx;
    if nvx && any(BayesInfo.p3(1:nvx)<0) 
        warning('Your prior allows for negative standard deviations for structural shocks. Due to working with variances, Dynare will be able to continue, but it is recommended to change your prior.')
    end
    offset=nvx;
    nvn=EstimatedParameters.nvn;
    if nvn && any(BayesInfo.p3(1+offset:offset+nvn)<0) 
        warning('Your prior allows for negative standard deviations for measurement error. Due to working with variances, Dynare will be able to continue, but it is recommended to change your prior.')
    end
    offset = nvx+nvn;
    ncx=EstimatedParameters.ncx; 
    if ncx && (any(BayesInfo.p3(1+offset:offset+ncx)<-1) || any(BayesInfo.p4(1+offset:offset+ncx)>1)) 
        warning('Your prior allows for correlations between structural shocks larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
    end
    offset = nvx+nvn+ncx;
    ncn=EstimatedParameters.ncn; 
    if ncn && (any(BayesInfo.p3(1+offset:offset+ncn)<-1) || any(BayesInfo.p4(1+offset:offset+ncn)>1)) 
        warning('Your prior allows for correlations between measurement errors larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
    end
end

% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(Model);

[lnprior, junk1,junk2,info]= priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
if info
    fprintf('The prior density evaluated at the initial values is Inf for the following parameters: %s\n',BayesInfo.name{info,1})
    error('The initial value of the prior is -Inf')
end

if DynareOptions.ramsey_policy
    %test whether specification matches
    inst_nbr = size(DynareOptions.instruments,1);
    if inst_nbr~=0
        orig_endo_aux_nbr = Model.orig_endo_nbr + min(find([Model.aux_vars.type] == 6)) - 1;
        implied_inst_nbr = orig_endo_aux_nbr - Model.orig_eq_nbr;
        if inst_nbr>implied_inst_nbr
            error('You have specified more instruments than there are omitted equations')
        elseif inst_nbr<implied_inst_nbr
            error('You have specified fewer instruments than there are omitted equations')
        end
    end
end

% Evaluate the likelihood.
ana_deriv = DynareOptions.analytic_derivation;
DynareOptions.analytic_derivation=0;
if ~isequal(DynareOptions.mode_compute,11) || ...
        (isequal(DynareOptions.mode_compute,11) && isequal(DynareOptions.order,1))
    %shut off potentially automatic switch to diffuse filter for the
    %purpose of checking stochastic singularity
    use_univariate_filters_if_singularity_is_detected_old=DynareOptions.use_univariate_filters_if_singularity_is_detected;
    DynareOptions.use_univariate_filters_if_singularity_is_detected=0;
    [fval,info] = feval(objective_function,xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults);
    if info(1)==50
        fprintf('\ninitial_estimation_checks:: The forecast error variance in the multivariate Kalman filter became singular.\n')
        fprintf('initial_estimation_checks:: This is often a sign of stochastic singularity, but can also sometimes happen by chance\n')
        fprintf('initial_estimation_checks:: for a particular combination of parameters and data realizations.\n')
        fprintf('initial_estimation_checks:: If you think the latter is the case, you should try with different initial values for the estimated parameters.\n')
        error('initial_estimation_checks:: The forecast error variance in the multivariate Kalman filter became singular.')
    end
    %reset options 
    DynareOptions.use_univariate_filters_if_singularity_is_detected=use_univariate_filters_if_singularity_is_detected_old;
else 
    info=0;
    fval = 0;
end
if DynareOptions.debug
    DynareResults.likelihood_at_initial_parameters=fval;
end
DynareOptions.analytic_derivation=ana_deriv;

% if DynareOptions.mode_compute==5
%     if ~strcmp(func2str(objective_function),'dsge_likelihood')
%         error('Options mode_compute=5 is not compatible with non linear filters or Dsge-VAR models!')
%     end
% end
if isnan(fval)
    error('The initial value of the likelihood is NaN')
elseif imag(fval)
    error('The initial value of the likelihood is complex')
end

if info(1) > 0
    disp('Error in computing likelihood for initial parameter values')
    print_info(info, DynareOptions.noprint, DynareOptions)
end

if DynareOptions.prefilter==1
    if (~DynareOptions.loglinear && any(abs(DynareResults.steady_state(BayesInfo.mfys))>1e-9)) || (DynareOptions.loglinear && any(abs(log(DynareResults.steady_state(BayesInfo.mfys)))>1e-9))
        disp(['You are trying to estimate a model with a non zero steady state for the observed endogenous'])
        disp(['variables using demeaned data!'])
        error('You should change something in your mod file...')
    end
end

if ~isequal(DynareOptions.mode_compute,11)
    disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
end
