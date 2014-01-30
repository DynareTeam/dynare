function [dataset_, xparam1, hh, M_, options_, oo_, estim_params_,bayestopt_, fake] = dynare_estimation_init(var_list_, dname, gsa_flag, M_, options_, oo_, estim_params_, bayestopt_)

% function dynare_estimation_init(var_list_, gsa_flag)
% preforms initialization tasks before estimation or
% global sensitivity analysis
%
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%   gsa_flag:   flag for GSA operation (optional)
%
% OUTPUTS
%   data:    data after required transformation
%   rawdata:  data as in the data file
%   xparam1:    initial value of estimated parameters as returned by
%               set_prior()
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2013 Dynare Team
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

hh = [];

if isempty(gsa_flag)
    gsa_flag = 0;
else% Decide if a DSGE or DSGE-VAR has to be estimated.
    if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
        options_.dsge_var = 1;
    end
    var_list_ = check_list_of_variables(options_, M_, var_list_);
    options_.varlist = var_list_;
end

% Get the indices of the observed variables in M_.endo_names.
options_.lgyidx2varobs = zeros(size(M_.endo_names,1),1);
for i = 1:size(M_.endo_names,1)
    tmp = strmatch(deblank(M_.endo_names(i,:)),options_.varobs,'exact');
    if ~isempty(tmp)
        if length(tmp)>1
            skipline()
            error(['Multiple declarations of ' deblank(M_.endo_names(i,:)) ' as an observed variable is not allowed!'])
        end
        options_.lgyidx2varobs(i) = tmp;
    end
end

if options_.order>2
    error(['I cannot estimate a model with a ' int2str(options_.order) ' order approximation of the model!'])
end

% Set options_.lik_init equal to 3 if diffuse filter is used or
% kalman_algo refers to a diffuse filter algorithm.
if (options_.diffuse_filter==1) || (options_.kalman_algo > 2)
    if options_.lik_init == 2
        error(['options diffuse_filter, lik_init and/or kalman_algo have ' ...
               'contradictory settings'])
    else
        options_.lik_init = 3;
    end
end

% If options_.lik_init == 1
%  set by default options_.qz_criterium to 1-1e-6
%  and check options_.qz_criterium < 1-eps if options_.lik_init == 1
% Else set by default options_.qz_criterium to 1+1e-6
if options_.lik_init == 1
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1-1e-6;
    elseif options_.qz_criterium > 1-eps
        error(['estimation: option qz_criterium is too large for estimating ' ...
               'a stationary model. If your model contains unit roots, use ' ...
               'option diffuse_filter'])
    end
elseif isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

% Set options related to filtered variables.
if ~isequal(options_.filtered_vars,0) && isempty(options_.filter_step_ahead)
    options_.filter_step_ahead = 1;
end
if ~isequal(options_.filtered_vars,0) && isequal(options_.filter_step_ahead,0)
    options_.filter_step_ahead = 1;
end
if ~isequal(options_.filter_step_ahead,0)
    options_.nk = max(options_.filter_step_ahead);
end

% Set the name of the directory where (intermediary) results will be saved.
if isempty(dname)
    M_.dname = M_.fname;
else
    M_.dname = dname;
end

% Set the number of observed variables.
n_varobs = size(options_.varobs,1);

% Set priors over the estimated parameters.
if ~isempty(estim_params_)
    [xparam1,estim_params_,bayestopt_,lb,ub,M_] = set_prior(estim_params_,M_,options_);
end

% Check that the provided mode_file is compatible with the current estimation settings.
if ~isempty(estim_params_) && ~isempty(options_.mode_file) && ~options_.mh_posterior_mode_estimation
    number_of_estimated_parameters = length(xparam1);
    mode_file = load(options_.mode_file);
    if number_of_estimated_parameters>length(mode_file.xparam1)
        % More estimated parameters than parameters in the mode file.
        skipline()
        disp(['The posterior mode file ' options_.mode_file ' has been generated using another specification of the model or another model!'])
        disp(['Your mode file contains estimates for ' int2str(length(mode_file.xparam1)) ' parameters, while you are attempting to estimate ' int2str(number_of_estimated_parameters) ' parameters:'])
        md = []; xd = [];
        for i=1:number_of_estimated_parameters
            id = strmatch(deblank(bayestopt_.name(i,:)),mode_file.parameter_names,'exact');
            if isempty(id)
                disp(['--> Estimated parameter ' bayestopt_.name{i} ' is not present in the loaded mode file (prior mean will be used, if possible).'])
            else
                xd = [xd; i];
                md = [md; id];
            end
        end
        for i=1:length(mode_file.xparam1)
            id = strmatch(mode_file.parameter_names{i},bayestopt_.name,'exact');
            if isempty(id)
                disp(['--> Parameter ' mode_file.parameter_names{i} ' is not estimated according to the current mod file.'])
            end
        end
        if ~options_.mode_compute
            % The posterior mode is not estimated.
            error('Please change the mode_file option, the list of estimated parameters or set mode_compute>0.')
        else
            % The posterior mode is estimated, the Hessian evaluated at the mode is not needed so we set values for the parameters missing in the mode file using the prior mean. 
            if ~isempty(xd)
                xparam1(xd) = mode_file.xparam1(md);
            else
                error('Please remove the mode_file option.')
            end
        end
    elseif number_of_estimated_parameters<length(mode_file.xparam1)
        % Less estimated parameters than parameters in the mode file.
        skipline()
        disp(['The posterior mode file ' options_.mode_file ' has been generated using another specification of the model or another model!'])
        disp(['Your mode file contains estimates for ' int2str(length(mode_file.xparam1)) ' parameters, while you are attempting to estimate only ' int2str(number_of_estimated_parameters) ' parameters:'])
        md = []; xd = [];
        for i=1:number_of_estimated_parameters
            id = strmatch(deblank(bayestopt_.name(i,:)),mode_file.parameter_names,'exact');
            if isempty(id)
                disp(['--> Estimated parameter ' deblank(bayestopt_.name(i,:)) ' is not present in the loaded mode file (prior mean will be used, if possible).'])
            else
                xd = [xd; i];
                md = [md; id];
            end
        end
        for i=1:length(mode_file.xparam1)
            id = strmatch(mode_file.parameter_names{i},bayestopt_.name,'exact');
            if isempty(id)
                disp(['--> Parameter ' mode_file.parameter_names{i} ' is not estimated according to the current mod file.'])
            end
        end
        if ~options_.mode_compute
            % The posterior mode is not estimated. If possible, fix the mode_file.
            if isequal(length(xd),number_of_estimated_parameters)
                disp('==> Fix mode file (remove unused parameters).')
                xparam1 = mode_file.xparam1(md);
                if isfield(mode_file,'hh')
                    hh = mode_file.hh(md,md);
                end
            else
                error('Please change the mode_file option, the list of estimated parameters or set mode_compute>0.')
            end
        else
            % The posterior mode is estimated, the Hessian evaluated at the mode is not needed so we set values for the parameters missing in the mode file using the prior mean. 
            if ~isempty(xd)
                xparam1(xd) = mode_file.xparam1(md);
            else
                % None of the estimated parameters are present in the mode_file.
                error('Please remove the mode_file option.')
            end
        end
    else
        % The number of declared estimated parameters match the number of parameters in the mode file. 
        % Check that the parameters in the mode file and according to the current mod file are identical.
        if isequal(mode_file.parameter_names, bayestopt_.name)
            xparam1 = mode_file.xparam1;
            if isfield(mode_file,'hh')
                hh = mode_file.hh;
            end
        else
            skipline()
            disp(['The posterior mode file ' options_.mode_file ' has been generated using another specification of the model or another model!'])
            % Check if this only an ordering issue or if the missing parameters can be initialized with the prior mean.
            md = []; xd = [];
            for i=1:number_of_estimated_parameters
                id = strmatch(deblank(bayestopt_.name(i,:)), mode_file.parameter_names,'exact');
                if isempty(id)
                    disp(['--> Estimated parameter ' bayestopt_.name{i} ' is not present in the loaded mode file.'])
                else
                    xd = [xd; i];
                    md = [md; id];
                end
            end
            if ~options_.mode_compute
                % The posterior mode is not estimated
                if isequal(length(xd), number_of_estimated_parameters)
                    % This is an ordering issue.
                    xparam1 = mode_file.xparam1(md);
                    if isfield(mode_file,'hh')
                        hh = mode_file.hh(md,md);
                    end
                else
                    error('Please change the mode_file option, the list of estimated parameters or set mode_compute>0.')
                end
            else
                % The posterior mode is estimated, the Hessian evaluated at the mode is not needed so we set values for the parameters missing in the mode file using the prior mean. 
                if ~isempty(xd)
                    xparam1(xd) = mode_file.xparam1(md);
                    if isfield(mode_file,'hh')
                        hh(xd,xd) = mode_file.hh(md,md);
                    end
                else
                    % None of the estimated parameters are present in the mode_file.
                    error('Please remove the mode_file option.')
                end
            end
        end
    end
    skipline()
end

if ~isempty(estim_params_) 
    if ~isempty(bayestopt_) && any(bayestopt_.pshape > 0)
        % Plot prior densities.
        if ~options_.nograph && options_.plot_priors
            plot_priors(bayestopt_,M_,estim_params_,options_)
        end
        % Set prior bounds
        bounds = prior_bounds(bayestopt_,options_);
        bounds(:,1)=max(bounds(:,1),lb);
        bounds(:,2)=min(bounds(:,2),ub);
    else  % estimated parameters but no declared priors
        % No priors are declared so Dynare will estimate the model by
        % maximum likelihood with inequality constraints for the parameters.
        options_.mh_replic = 0;% No metropolis.
        bounds(:,1) = lb;
        bounds(:,2) = ub;
    end
    % Test if initial values of the estimated parameters are all between the prior lower and upper bounds.
    check_prior_bounds(xparam1,bounds,M_,estim_params_,options_,bayestopt_)

    lb = bounds(:,1);
    ub = bounds(:,2);
    bayestopt_.lb = lb;
    bayestopt_.ub = ub;
end

if isempty(estim_params_)% If estim_params_ is empty (e.g. when running the smoother on a calibrated model)
    if ~options_.smoother
        error('ESTIMATION: the ''estimated_params'' block is mandatory (unless you are running a smoother)')
    end
    xparam1 = [];
    bayestopt_.lb = [];
    bayestopt_.ub = [];
    bayestopt_.jscale = [];
    bayestopt_.pshape = [];
    bayestopt_.p1 = [];
    bayestopt_.p2 = [];
    bayestopt_.p3 = [];
    bayestopt_.p4 = [];
    bayestopt_.p5 = [];
    bayestopt_.p6 = [];
    bayestopt_.p7 = [];
    estim_params_.nvx = 0;
    estim_params_.nvn = 0;
    estim_params_.ncx = 0;
    estim_params_.ncn = 0;
    estim_params_.np = 0;
end

% storing prior parameters in results
oo_.prior.mean = bayestopt_.p1;
oo_.prior.variance = diag(bayestopt_.p2.^2);

% Is there a linear trend in the measurement equation?
if ~isfield(options_,'trend_coeffs') % No!
    bayestopt_.with_trend = 0;
else% Yes!
    bayestopt_.with_trend = 1;
    bayestopt_.trend_coeff = {};
    trend_coeffs = options_.trend_coeffs;
    nt = length(trend_coeffs);
    for i=1:n_varobs
        if i > length(trend_coeffs)
            bayestopt_.trend_coeff{i} = '0';
        else
            bayestopt_.trend_coeff{i} = trend_coeffs{i};
        end
    end
end

% Set the "size" of penalty.
objective_function_penalty_base = 1e8;

% Get informations about the variables of the model.
dr = set_state_space(oo_.dr,M_,options_);
oo_.dr = dr;
nstatic = M_.nstatic;          % Number of static variables.
npred = M_.nspred;             % Number of predetermined variables.
nspred = M_.nspred;            % Number of predetermined variables in the state equation.

% Test if observed variables are declared.
if isempty(options_.varobs)
    error('VAROBS is missing')
end

% Setting resticted state space (observed + predetermined variables)
var_obs_index = [];
k1 = [];
for i=1:n_varobs
    var_obs_index = [var_obs_index; strmatch(deblank(options_.varobs(i,:)),M_.endo_names(dr.order_var,:),'exact')];
    k1 = [k1; strmatch(deblank(options_.varobs(i,:)),M_.endo_names, 'exact')];
end

% Define union of observed and state variables
k2 = union(var_obs_index,[M_.nstatic+1:M_.nstatic+M_.nspred]', 'rows');
% Set restrict_state to postion of observed + state variables in expanded state vector.
oo_.dr.restrict_var_list = k2;
bayestopt_.restrict_var_list = k2;
% set mf0 to positions of state variables in restricted state vector for likelihood computation.
[junk,bayestopt_.mf0] = ismember([M_.nstatic+1:M_.nstatic+M_.nspred]',k2);
% Set mf1 to positions of observed variables in restricted state vector for likelihood computation.
[junk,bayestopt_.mf1] = ismember(var_obs_index,k2);
% Set mf2 to positions of observed variables in expanded state vector for filtering and smoothing.
bayestopt_.mf2  = var_obs_index;
bayestopt_.mfys = k1;

[junk,ic] = intersect(k2,nstatic+(1:npred)');
oo_.dr.restrict_columns = [ic; length(k2)+(1:nspred-npred)'];

k3 = [];
k3p = [];
if options_.selected_variables_only
    for i=1:size(var_list_,1)
        k3 = [k3; strmatch(var_list_(i,:),M_.endo_names(dr.order_var,:), ...
                           'exact')];
        k3p = [k3; strmatch(var_list_(i,:),M_.endo_names, ...
                           'exact')];
    end
else
    k3 = (1:M_.endo_nbr)';
    k3p = (1:M_.endo_nbr)';
end

% Define union of observed and state variables
if options_.block == 1
    k1 = k1';
    [k2, i_posA, i_posB] = union(k1', M_.state_var', 'rows');
    % Set restrict_state to postion of observed + state variables in expanded state vector.
    oo_.dr.restrict_var_list  = [k1(i_posA) M_.state_var(sort(i_posB))];
    % set mf0 to positions of state variables in restricted state vector for likelihood computation.
    [junk,bayestopt_.mf0] = ismember(M_.state_var',oo_.dr.restrict_var_list);
    % Set mf1 to positions of observed variables in restricted state vector for likelihood computation.
    [junk,bayestopt_.mf1] = ismember(k1,oo_.dr.restrict_var_list);
    % Set mf2 to positions of observed variables in expanded state vector for filtering and smoothing.
    bayestopt_.mf2  = var_obs_index;
    bayestopt_.mfys = k1;
    oo_.dr.restrict_columns = [size(i_posA,1)+(1:size(M_.state_var,2))];

    [k2, i_posA, i_posB] = union(k3p, M_.state_var', 'rows');
    bayestopt_.smoother_var_list = [k3p(i_posA); M_.state_var(sort(i_posB))'];
    [junk,junk,bayestopt_.smoother_saved_var_list] = intersect(k3p,bayestopt_.smoother_var_list(:));
    [junk,ic] = intersect(bayestopt_.smoother_var_list,M_.state_var);
    bayestopt_.smoother_restrict_columns = ic;
    [junk,bayestopt_.smoother_mf] = ismember(k1, ...
                                             bayestopt_.smoother_var_list);
else
    k2 = union(var_obs_index,[M_.nstatic+1:M_.nstatic+M_.nspred]', 'rows');
    % Set restrict_state to postion of observed + state variables in expanded state vector.
    oo_.dr.restrict_var_list = k2;
    % set mf0 to positions of state variables in restricted state vector for likelihood computation.
    [junk,bayestopt_.mf0] = ismember([M_.nstatic+1:M_.nstatic+M_.nspred]',k2);
    % Set mf1 to positions of observed variables in restricted state vector for likelihood computation.
    [junk,bayestopt_.mf1] = ismember(var_obs_index,k2);
    % Set mf2 to positions of observed variables in expanded state vector for filtering and smoothing.
    bayestopt_.mf2  = var_obs_index;
    bayestopt_.mfys = k1;
    [junk,ic] = intersect(k2,nstatic+(1:npred)');
    oo_.dr.restrict_columns = [ic; length(k2)+(1:nspred-npred)'];

    bayestopt_.smoother_var_list = union(k2,k3);
    [junk,junk,bayestopt_.smoother_saved_var_list] = intersect(k3,bayestopt_.smoother_var_list(:));
    [junk,ic] = intersect(bayestopt_.smoother_var_list,nstatic+(1:npred)');
    bayestopt_.smoother_restrict_columns = ic;
    [junk,bayestopt_.smoother_mf] = ismember(var_obs_index, ...
                                             bayestopt_.smoother_var_list);
end;

if options_.analytic_derivation,
    options_.analytic_derivation = 1;
    if ~(exist('sylvester3','file')==2),
        dynareroot = strrep(which('dynare'),'dynare.m','');
        addpath([dynareroot 'gensylv'])
    end
    if estim_params_.np,
        % check if steady state changes param values
        M=M_;
        M.params(estim_params_.param_vals(:,1)) = M.params(estim_params_.param_vals(:,1))*1.01;
        if options_.diffuse_filter
            steadystate_check_flag = 0;
        else
            steadystate_check_flag = 1;
        end
        [tmp1, params] = evaluate_steady_state(oo_.steady_state,M,options_,oo_,steadystate_check_flag);
        change_flag=any(find(params-M.params));
        if change_flag,
            disp('The steadystate file changed the values for the following parameters: '),
            disp(M.param_names(find(params-M.params),:))
            disp('The derivatives of jacobian and steady-state will be computed numerically'),
            disp('(re-set options_.analytic_derivation_mode= -2)'),
            options_.analytic_derivation_mode= -2;
        end
    end
end

% Test if the data file is declared.
if isempty(options_.datafile)
    if gsa_flag
        dataset_ = [];
%         rawdata = [];
%         data_info = [];
        return
    else
        error('datafile option is missing')
    end
end

% If jscale isn't specified for an estimated parameter, use global option options_.jscale, set to 0.2, by default.
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

% Load and transform data.
transformation = [];
if options_.loglinear && ~options_.logdata
    transformation = @log;
end
xls.sheet = options_.xls_sheet;
xls.range = options_.xls_range;

if ~isfield(options_,'nobs')
    options_.nobs = [];
end

dataset_ = initialize_dataset(options_.datafile,options_.varobs,options_.first_obs,options_.nobs,transformation,options_.prefilter,xls);

options_.nobs = dataset_.info.ntobs;

% setting steadystate_check_flag option
if options_.diffuse_filter
    steadystate_check_flag = 0;
else
    steadystate_check_flag = 1;
end

M = M_;
nvx = estim_params_.nvx;
ncx = estim_params_.ncx;
nvn = estim_params_.nvn;
ncn = estim_params_.ncn;
if estim_params_.np
  M.params(estim_params_.param_vals(:,1)) = xparam1(nvx+ncx+nvn+ncn+1:end);
end
[oo_.steady_state, params] = evaluate_steady_state(oo_.steady_state,M,options_,oo_,steadystate_check_flag);
if all(abs(oo_.steady_state(bayestopt_.mfys))<1e-9)
    options_.noconstant = 1;
else
    options_.noconstant = 0;
    % If the data are prefiltered then there must not be constants in the
    % measurement equation of the DSGE model or in the DSGE-VAR model.
    if options_.prefilter == 1
        fprintf('\nestimation_init: You have specified the option "prefilter" to demean your data,\n')
        fprintf('estimation_init: but your observation equations are not mean zero. Either change your observation\n')
        fprintf('estimation_init: equation or drop the prefiltering.\n')
        error('The option "prefilter" is inconsistent with the non-zero mean measurement equations in the model.')
    end
end

