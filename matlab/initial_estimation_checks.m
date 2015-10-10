function DynareResults = initial_estimation_checks(objective_function,xparam1,DynareDataset,DatasetInfo,Model,EstimatedParameters,DynareOptions,BayesInfo,BoundsInfo,DynareResults)
% function initial_estimation_checks(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)
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
%   DynareResults       [structure] storing the results
%
% OUTPUTS
%    DynareResults     structure of temporary results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2015 Dynare Team
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

if DynareDataset.vobs>Model.exo_nbr+EstimatedParameters.nvn
    error(['initial_estimation_checks:: Estimation can''t take place because there are less declared shocks than observed variables!'])
end

if DynareDataset.vobs>length(find(diag(Model.Sigma_e)))+EstimatedParameters.nvn
    error(['initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance!'])
end

if ~(strcmpi(DynareOptions.proposal_distribution, 'rand_multivariate_student') || ...
     strcmpi(DynareOptions.proposal_distribution, 'rand_multivariate_normal'))
    error(['initial_estimation_checks:: the proposal_distribution option to estimation takes either ' ...
        'rand_multivariate_student or rand_multivariate_normal as options']);
end

if DynareOptions.student_degrees_of_freedom <= 0
    error('initial_estimation_checks:: the student_degrees_of_freedom takes a positive integer argument');
end

if DynareOptions.TaRB.use_TaRB && (DynareOptions.TaRB.new_block_probability<0 || DynareOptions.TaRB.new_block_probability>1)
    error(['initial_estimation_checks:: The tarb_new_block_probability must be between 0 and 1!'])
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
        disp(Model.param_names(changed_par_indices,:));
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
        warning('Your prior allows for correlations between measurement errors larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
    end
    offset = nvx+nvn+ncx;
    ncn=EstimatedParameters.ncn; 
    if ncn && (any(BayesInfo.p3(1+offset:offset+ncn)<-1) || any(BayesInfo.p4(1+offset:offset+ncn)>1)) 
        warning('Your prior allows for correlations between structural shocks larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
    end
end

% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(Model);

[lnprior, junk1,junk2,info]= priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
if info
    fprintf('The prior density evaluated at the initial values is Inf for the following parameters: %s\n',BayesInfo.name{info,1})
    error('The initial value of the prior is -Inf')
end
    
% Evaluate the likelihood.
ana_deriv = DynareOptions.analytic_derivation;
DynareOptions.analytic_derivation=0;
if ~isequal(DynareOptions.mode_compute,11) || ...
        (isequal(DynareOptions.mode_compute,11) && isequal(DynareOptions.order,1))
  [fval,junk1,junk2,a,b,c,d] = feval(objective_function,xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults);
else 
    b=0;
    fval = 0;
end
if DynareOptions.debug
    DynareResults.likelihood_at_initial_parameters=fval;
end
DynareOptions.analytic_derivation=ana_deriv;

if DynareOptions.dsge_var || strcmp(func2str(objective_function),'non_linear_dsge_likelihood')
    info = b;
else
    info = d;
end

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
