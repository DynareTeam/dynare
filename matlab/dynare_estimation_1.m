function dynare_estimation_1(var_list_,dname)
% function dynare_estimation_1(var_list_,dname)
% runs the estimation of the model
%
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info

% Set particle filter flag.

if ~isfield( options_, 'non_central_approximation' )
    options_.non_central_approximation = 0;
end
if ~isfield( options_, 'gaussian_approximation' )
    options_.gaussian_approximation = 0;
end
if options_.order > 1 && ( options_.non_central_approximation == 0 ) && ( options_.gaussian_approximation == 0 )
    if options_.particle.status && options_.order==2
        skipline()
        disp('Estimation using a non linear filter!')
        skipline()
        if ~options_.nointeractive && ismember(options_.mode_compute,[1,3,4]) && ~strcmpi(options_.particle.filter_algorithm,'gf')% Known gradient-based optimizers
            disp('You are using a gradient-based mode-finder. Particle filtering introduces discontinuities in the') 
            disp('objective function w.r.t the parameters. Thus, should use a non-gradient based optimizer.')
            fprintf('\nPlease choose a mode-finder:\n')
            fprintf('\t 0 - Continue using gradient-based method (it is most likely that you will no get any sensible result).\n')
            fprintf('\t 6 - Monte Carlo based algorithm\n')
            fprintf('\t 7 - Nelder-Mead simplex based optimization routine (Matlab optimization toolbox required)\n')
            fprintf('\t 8 - Nelder-Mead simplex based optimization routine (Dynare''s implementation)\n')
            fprintf('\t 9 - CMA-ES (Covariance Matrix Adaptation Evolution Strategy) algorithm\n')
            choice = [];
            while isempty(choice)
                choice = input('Please enter your choice: ');
                if isnumeric(choice) && isint(choice) && ismember(choice,[0 6 7 8 9])
                    if choice
                        options_.mode_compute = choice;
                    end
                else
                    fprintf('\nThis is an invalid choice (you have to choose between 0, 6, 7, 8 and 9).\n')
                    choice = [];
                end
            end
        end
    elseif options_.particle.status && options_.order>2
        error(['Non linear filter are not implemented with order ' int2str(options_.order) ' approximation of the model!'])
    elseif ~options_.particle.status && options_.order==2
        error('For estimating the model with a second order approximation using a non linear filter, one should have options_.particle.status=1;')
    else
        error(['Cannot estimate a model with an order ' int2str(options_.order) ' approximation!'])
    end
end

if ~options_.dsge_var
    if options_.particle.status && ~( options_.non_central_approximation || options_.gaussian_approximation )
        objective_function = str2func('non_linear_dsge_likelihood');
        if strcmpi(options_.particle.filter_algorithm, 'sis')
            options_.particle.algorithm = 'sequential_importance_particle_filter';
        elseif strcmpi(options_.particle.filter_algorithm, 'apf')
            options_.particle.algorithm = 'auxiliary_particle_filter';
        elseif strcmpi(options_.particle.filter_algorithm, 'gf')
            options_.particle.algorithm = 'gaussian_filter';
        elseif strcmpi(options_.particle.filter_algorithm,  'gmf')
            options_.particle.algorithm = 'gaussian_mixture_filter';
        elseif strcmpi(options_.particle.filter_algorithm, 'cpf')
            options_.particle.algorithm = 'conditional_particle_filter';
        else
            error(['Estimation: Unknown filter ' options_.particle.filter_algorithm])
        end
    else
        objective_function = str2func('dsge_likelihood');
    end
else
    objective_function = str2func('dsge_var_likelihood');
end

[dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_, bayestopt_, bounds] = ...
    dynare_estimation_init(var_list_, dname, [], M_, options_, oo_, estim_params_, bayestopt_);

if options_.dsge_var
    check_dsge_var_model(M_, estim_params_, bayestopt_);
    if dataset_info.missing.state
        error('Estimation::DsgeVarLikelihood: I cannot estimate a DSGE-VAR model with missing observations!')
    end
    if options_.noconstant
        var_sample_moments(options_.dsge_varlag, -1, dataset_);
    else
        % The steady state is non zero ==> a constant in the VAR is needed!
        var_sample_moments(options_.dsge_varlag, 0, dataset_);
    end
end

% Set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file).
M_.sigma_e_is_diagonal = 1;
if estim_params_.ncx || any(nnz(tril(M_.Correlation_matrix,-1))) || isfield(estim_params_,'calibrated_covariances')
    M_.sigma_e_is_diagonal = 0;
end

data = dataset_.data;
rawdata = dataset_info.rawdata;
data_index = dataset_info.missing.aindex;
missing_value = dataset_info.missing.state;

% Set number of observations
gend = dataset_.nobs;
% Set the number of observed variables.
n_varobs = length(options_.varobs);
% Get the number of parameters to be estimated.
nvx = estim_params_.nvx;  % Variance of the structural innovations (number of parameters).
nvn = estim_params_.nvn;  % Variance of the measurement innovations (number of parameters).
ncx = estim_params_.ncx;  % Covariance of the structural innovations (number of parameters).
ncn = estim_params_.ncn;  % Covariance of the measurement innovations (number of parameters).
np  = estim_params_.np ;  % Number of deep parameters.
nx  = nvx+nvn+ncx+ncn+np; % Total number of parameters to be estimated.
%% Set the names of the priors.
pnames = ['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];

dr = oo_.dr;

% load optimal_mh_scale parameter if previous run was with mode_compute=6
mh_scale_fname = [M_.fname '_optimal_mh_scale_parameter.mat'];
if exist(mh_scale_fname)
    if options_.mode_compute == 0
        tmp = load(mh_scale_fname,'Scale');
        bayestopt_.mh_jscale = tmp.Scale;
        clear tmp;
    else
        % remove the file if mode_compute ~= 0
        delete(mh_scale_fname)
    end
end

if ~isempty(estim_params_)
    M_ = set_all_parameters(xparam1,estim_params_,M_);
end


%% perform initial estimation checks;
try
    oo_ = initial_estimation_checks(objective_function,xparam1,dataset_,dataset_info,M_,estim_params_,options_,bayestopt_,bounds,oo_);
catch % if check fails, provide info on using calibration if present
    e = lasterror();
    if estim_params_.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('ESTIMATION_CHECKS: There was an error in computing the likelihood for initial parameter values.\n')
        fprintf('ESTIMATION_CHECKS: You should try using the calibrated version of the model as starting values. To do\n')
        fprintf('ESTIMATION_CHECKS: this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n')    
        fprintf('ESTIMATION_CHECKS: command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(e);
end

if isequal(options_.mode_compute,0) && isempty(options_.mode_file) && options_.mh_posterior_mode_estimation==0
    if options_.smoother == 1
        [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,gend,transpose(data),data_index,missing_value);
        oo_.Smoother.SteadyState = ys;
        oo_.Smoother.TrendCoeffs = trend_coeff;
        if options_.filter_covariance
            oo_.Smoother.Variance = P;
        end
        i_endo = bayestopt_.smoother_saved_var_list;
        if options_.nk ~= 0
            oo_.FilteredVariablesKStepAhead = ...
                aK(options_.filter_step_ahead,i_endo,:);
            if ~isempty(PK)
                oo_.FilteredVariablesKStepAheadVariances = ...
                    PK(options_.filter_step_ahead,i_endo,i_endo,:);
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition = ...
                    decomp(options_.filter_step_ahead,i_endo,:,:);
            end
        end
        for i=bayestopt_.smoother_saved_var_list'
            i1 = dr.order_var(bayestopt_.smoother_var_list(i));
            eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ...
                  ' = atT(i,:)'';']);
            if options_.nk > 0
                eval(['oo_.FilteredVariables.' deblank(M_.endo_names(i1,:)) ...
                      ' = squeeze(aK(1,i,2:end-(options_.nk-1)));']);
            end
            eval(['oo_.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ' = updated_variables(i,:)'';']);
        end
        for i=1:M_.exo_nbr
            eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
        end
    end
    return
end

% Estimation of the posterior mode or likelihood mode
if ~isequal(options_.mode_compute,0) && ~options_.mh_posterior_mode_estimation
    %prepare settings for newrat
    if options_.mode_compute==5
        %get whether outer product Hessian is requested
        newratflag=[];
        if ~isempty(options_.optim_opt)
            options_list = read_key_value_string(options_.optim_opt);
            for i=1:rows(options_list)
                if strcmp(options_list{i,1},'Hessian')
                    newratflag=options_list{i,2};
                end
            end
        end
        if options_.analytic_derivation,
            options_analytic_derivation_old = options_.analytic_derivation;
            options_.analytic_derivation = -1;
            if ~isempty(newratflag) && newratflag~=0 %numerical hessian explicitly specified
                error('newrat: analytic_derivation is incompatible with numerical Hessian.')
            else %use default
                newratflag=0; %exclude DYNARE numerical hessian
            end
        elseif ~options_.analytic_derivation 
            if isempty(newratflag) 
                newratflag=options_.newrat.hess; %use default numerical dynare hessian                
            end
        end
    end
    
    [xparam1, fval, exitflag, hh, options_, Scale] = dynare_minimize_objective(objective_function,xparam1,options_.mode_compute,options_,[bounds.lb bounds.ub],bayestopt_.name,bayestopt_,hh,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
    fprintf('\nFinal value of minus the log posterior (or likelihood):%f \n', fval);

    if isnumeric(options_.mode_compute) && options_.mode_compute==5 && options_.analytic_derivation==-1 %reset options changed by newrat
            options_.analytic_derivation = options_analytic_derivation_old; %reset      
    elseif isnumeric(options_.mode_compute) && options_.mode_compute==6 %save scaling factor
        save([M_.fname '_optimal_mh_scale_parameter.mat'],'Scale');
        options_.mh_jscale = Scale;
        bayestopt_.jscale = ones(length(xparam1),1)*Scale;
    end       
    if ~isnumeric(options_.mode_compute) || ~isequal(options_.mode_compute,6) %always already computes covariance matrix
        if options_.cova_compute == 1 %user did not request covariance not to be computed
            if options_.analytic_derivation && strcmp(func2str(objective_function),'dsge_likelihood'),
                ana_deriv_old = options_.analytic_derivation;
                options_.analytic_derivation = 2;
                [junk1, junk2, hh] = feval(objective_function,xparam1, ...
                    dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
                options_.analytic_derivation = ana_deriv_old;
            elseif ~isnumeric(options_.mode_compute) || ~(isequal(options_.mode_compute,5) && newratflag~=1), 
                % with flag==0, we force to use the hessian from outer
                % product gradient of optimizer 5
                hh = reshape(hessian(objective_function,xparam1, ...
                    options_.gstep,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_),nx,nx);
            elseif isnumeric(options_.mode_compute) && isequal(options_.mode_compute,5)
                % other numerical hessian options available with optimizer 5
                %
                % if newratflag == 0
                % compute outer product gradient of optimizer 5
                %
                % if newratflag == 2
                % compute 'mixed' outer product gradient of optimizer 5
                % with diagonal elements computed with numerical second order derivatives
                %
                % uses univariate filters, so to get max # of available
                % densitities for outer product gradient
                kalman_algo0 = options_.kalman_algo;
                compute_hessian = 1;
                if ~((options_.kalman_algo == 2) || (options_.kalman_algo == 4)),
                    options_.kalman_algo=2;
                    if options_.lik_init == 3,
                        options_.kalman_algo=4;
                    end
                elseif newratflag==0, % hh already contains outer product gradient with univariate filter
                    compute_hessian = 0;                                            
                end
                if compute_hessian,
                    crit = options_.newrat.tolerance.f;
                    newratflag = newratflag>0;
                    hh = reshape(mr_hessian(0,xparam1,objective_function,newratflag,crit,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_), nx, nx);
                end
                options_.kalman_algo = kalman_algo0;
            end
        end
    end
    parameter_names = bayestopt_.name;
    if options_.cova_compute || options_.mode_compute==5 || options_.mode_compute==6
        save([M_.fname '_mode.mat'],'xparam1','hh','parameter_names','fval');
    else
        save([M_.fname '_mode.mat'],'xparam1','parameter_names','fval');
    end
end

switch options_.MCMC_jumping_covariance
    case 'hessian' %Baseline
        %do nothing and use hessian from mode_compute
    case 'prior_variance' %Use prior variance
        if any(isinf(bayestopt_.p2))
            error('Infinite prior variances detected. You cannot use the prior variances as the proposal density, if some variances are Inf.')
        else
            hh = diag(1./(bayestopt_.p2.^2));
        end
    case 'identity_matrix' %Use identity
        hh = eye(nx);
    otherwise %user specified matrix in file
        try
            load(options_.MCMC_jumping_covariance,'jumping_covariance')
            hh=jumping_covariance;
        catch
            error(['No matrix named ''jumping_covariance'' could be found in ',options_.MCMC_jumping_covariance,'.mat'])
        end
        [nrow, ncol]=size(hh);
        if ~isequal(nrow,ncol) && ~isequal(nrow,nx) %check if square and right size
            error(['jumping_covariance matrix must be square and have ',num2str(nx),' rows and columns'])
        end
        try %check for positive definiteness
            chol(hh);
        catch
            error(['Specified jumping_covariance is not positive definite'])
        end
end

if ~options_.mh_posterior_mode_estimation && options_.cova_compute
    try
        chol(hh);
    catch
        skipline()
        disp('POSTERIOR KERNEL OPTIMIZATION PROBLEM!')
        disp(' (minus) the hessian matrix at the "mode" is not positive definite!')
        disp('=> posterior variance of the estimated parameters are not positive.')
        disp('You should try to change the initial values of the parameters using')
        disp('the estimated_params_init block, or use another optimization routine.')
        params_at_bound=find(xparam1==bounds.ub | xparam1==bounds.lb);
        if ~isempty(params_at_bound)
            for ii=1:length(params_at_bound)
            params_at_bound_name{ii,1}=get_the_name(params_at_bound(ii),0,M_,estim_params_,options_);
            end
            disp_string=[params_at_bound_name{1,:}];
            for ii=2:size(params_at_bound_name,1)
                disp_string=[disp_string,', ',params_at_bound_name{ii,:}];
            end
            fprintf('\nThe following parameters are at the prior bound: %s\n', disp_string)
            fprintf('Some potential solutions are:\n')
            fprintf('   - Check your model for mistakes.\n')
            fprintf('   - Check whether model and data are consistent (correct observation equation).\n')
            fprintf('   - Shut off prior_trunc.\n')
            fprintf('   - Change the optimization bounds.\n')
            fprintf('   - Use a different mode_compute like 6 or 9.\n')
            fprintf('   - Check whether the parameters estimated are identified.\n')
            fprintf('   - Increase the informativeness of the prior.\n')
        end
        warning('The results below are most likely wrong!');
    end
end

if options_.mode_check.status == 1 && ~options_.mh_posterior_mode_estimation
    ana_deriv_old = options_.analytic_derivation;
    options_.analytic_derivation = 0;
    mode_check(objective_function,xparam1,hh,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
    options_.analytic_derivation = ana_deriv_old;
end

oo_.posterior.optimization.mode = xparam1;
oo_.posterior.optimization.Variance = [];
if ~options_.mh_posterior_mode_estimation
    if options_.cova_compute
        invhess = inv(hh);
        stdh = sqrt(diag(invhess));
        oo_.posterior.optimization.Variance = invhess;
    end
else
    variances = bayestopt_.p2.*bayestopt_.p2;
    idInf = isinf(variances);
    variances(idInf) = 1;
    invhess = options_.mh_posterior_mode_estimation*diag(variances);
    xparam1 = bayestopt_.p5;
    idNaN = isnan(xparam1);
    xparam1(idNaN) = bayestopt_.p1(idNaN);
end

if ~options_.cova_compute
    stdh = NaN(length(xparam1),1);
end

if any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    % display results table and store parameter estimates and standard errors in results
    oo_=display_estimation_results_table(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_,pnames,'Posterior','posterior');
    % Laplace approximation to the marginal log density:
    if options_.cova_compute
        estim_params_nbr = size(xparam1,1);
        scale_factor = -sum(log10(diag(invhess)));
        log_det_invhess = -estim_params_nbr*log(scale_factor)+log(det(scale_factor*invhess));
        likelihood = feval(objective_function,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
        oo_.MarginalDensity.LaplaceApproximation = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess - likelihood;
        skipline()
        disp(sprintf('Log data density [Laplace approximation] is %f.',oo_.MarginalDensity.LaplaceApproximation))
        skipline()
    end
elseif ~any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    oo_=display_estimation_results_table(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_,pnames,'Maximum Likelihood','mle');
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.fname '_params.mat'],'pindx');
end

if (any(bayestopt_.pshape  >0 ) && options_.mh_replic) || ...
        (any(bayestopt_.pshape >0 ) && options_.load_mh_file)  %% not ML estimation
    bounds = prior_bounds(bayestopt_,options_);
    outside_bound_pars=find(xparam1 < bounds.lb | xparam1 > bounds.ub);
    if ~isempty(outside_bound_pars)
        for ii=1:length(outside_bound_pars)
            outside_bound_par_names{ii,1}=get_the_name(ii,0,M_,estim_params_,options_);
        end
        disp_string=[outside_bound_par_names{1,:}];
        for ii=2:size(outside_bound_par_names,1)
            disp_string=[disp_string,', ',outside_bound_par_names{ii,:}];
        end
        if options_.prior_trunc>0
            error(['Estimation:: Mode value(s) of ', disp_string ,' are outside parameter bounds. Potentially, you should set prior_trunc=0.'])
        else
            error(['Estimation:: Mode value(s) of ', disp_string ,' are outside parameter bounds.'])
        end
    end
    % runs MCMC
    if options_.mh_replic
        if options_.load_mh_file && options_.use_mh_covariance_matrix
            invhess = compute_mh_covariance_matrix;
        end
        ana_deriv_old = options_.analytic_derivation;
        options_.analytic_derivation = 0;
        if options_.cova_compute
            feval(options_.posterior_sampling_method,objective_function,options_.proposal_distribution,xparam1,invhess,bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_);
        else
            error('I Cannot start the MCMC because the Hessian of the posterior kernel at the mode was not computed.')
        end
        options_.analytic_derivation = ana_deriv_old;
    end
    if options_.mh_posterior_mode_estimation
        CutSample(M_, options_, estim_params_);
        return
    else
        if ~options_.nodiagnostic && options_.mh_replic>0
            oo_= McMCDiagnostics(options_, estim_params_, M_,oo_);
        end
        %% Here i discard first half of the draws:
        CutSample(M_, options_, estim_params_);
        %% Estimation of the marginal density from the Mh draws:
        if options_.mh_replic
            [marginal,oo_] = marginal_density(M_, options_, estim_params_, oo_);
            % Store posterior statistics by parameter name
            oo_ = GetPosteriorParametersStatistics(estim_params_, M_, options_, bayestopt_, oo_);
            if ~options_.nograph
                oo_ = PlotPosteriorDistributions(estim_params_, M_, options_, bayestopt_, oo_);
            end
            % Store posterior mean in a vector and posterior variance in
            % a matrix
            [oo_.posterior.metropolis.mean,oo_.posterior.metropolis.Variance] ...
                = GetPosteriorMeanVariance(M_,options_.mh_drop);
        else
            load([M_.fname '_results'],'oo_');
        end
        error_flag = metropolis_draw(1);
        if options_.bayesian_irf
            if error_flag
                error('Estimation::mcmc: I cannot compute the posterior IRFs!')
            end
            PosteriorIRF('posterior');
        end
        if options_.moments_varendo
            if error_flag
                error('Estimation::mcmc: I cannot compute the posterior moments for the endogenous variables!')
            end
            oo_ = compute_moments_varendo('posterior',options_,M_,oo_,var_list_);
        end
        if options_.smoother || ~isempty(options_.filter_step_ahead) || options_.forecast
            if error_flag
                error('Estimation::mcmc: I cannot compute the posterior statistics!')
            end
            prior_posterior_statistics('posterior',dataset_,dataset_info);
        end
        xparam1 = get_posterior_parameters('mean');
        M_ = set_all_parameters(xparam1,estim_params_,M_);
    end
end

if options_.particle.status
    return
end

if (~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape ...
                                                      > 0) && options_.load_mh_file)) ...
    || ~options_.smoother ) && ( options_.partial_information == 0 ) && ( options_.gaussian_approximation == 0 )% to be fixed
    %% ML estimation, or posterior mode without metropolis-hastings or metropolis without bayesian smooth variable
    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,dataset_.nobs,transpose(dataset_.data),dataset_info.missing.aindex,dataset_info.missing.state);
    oo_.Smoother.SteadyState = ys;
    oo_.Smoother.TrendCoeffs = trend_coeff;
    oo_.Smoother.Variance = P;
    i_endo = bayestopt_.smoother_saved_var_list;
    if ~isempty(options_.nk) && options_.nk ~= 0 && (~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file)))
        oo_.FilteredVariablesKStepAhead = aK(options_.filter_step_ahead, ...
                                             i_endo,:);
        if isfield(options_,'kalman_algo')
            if ~isempty(PK)
                oo_.FilteredVariablesKStepAheadVariances = ...
                    PK(options_.filter_step_ahead,i_endo,i_endo,:);
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition = ...
                    decomp(options_.filter_step_ahead,i_endo,:,:);
            end
        end
    end
    for i=bayestopt_.smoother_saved_var_list'
        i1 = dr.order_var(bayestopt_.smoother_var_list(i));
        eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ' = ' ...
                            'atT(i,:)'';']);
        if ~isempty(options_.nk) && options_.nk > 0 && ~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file))
            eval(['oo_.FilteredVariables.' deblank(M_.endo_names(i1,:)) ...
                  ' = squeeze(aK(1,i,2:end-(options_.nk-1)));']);
        end
        eval(['oo_.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ...
              ' = updated_variables(i,:)'';']);
    end
    for i=1:M_.exo_nbr
        eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
    end
    if ~options_.nograph,
        [nbplt,nr,nc,lr,lc,nstar] = pltorg(M_.exo_nbr);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([M_.fname '_SmoothedShocks.TeX'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation_1.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        for plt = 1:nbplt,
            fh = dyn_figure(options_,'Name','Smoothed shocks');
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            nstar0=min(nstar,M_.exo_nbr-(plt-1)*nstar);
            if gend==1
                marker_string{1,1}='-ro';
                marker_string{2,1}='-ko';
            else
                marker_string{1,1}='-r';
                marker_string{2,1}='-k';
            end
            for i=1:nstar0,
                k = (plt-1)*nstar+i;
                subplot(nr,nc,i);
                plot([1 gend],[0 0],marker_string{1,1},'linewidth',.5)
                hold on
                plot(1:gend,innov(k,:),marker_string{2,1},'linewidth',1)
                hold off
                name = deblank(M_.exo_names(k,:));
                if isempty(NAMES)
                    NAMES = name;
                else
                    NAMES = char(NAMES,name);
                end
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                if gend>1
                    xlim([1 gend])
                end
                if options_.TeX
                    texname = M_.exo_names_tex(k,:);
                    if isempty(TeXNAMES)
                        TeXNAMES = ['$ ' deblank(texname) ' $'];
                    else
                        TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                    end
                end
                title(name,'Interpreter','none')
            end
            dyn_saveas(fh,[M_.fname '_SmoothedShocks' int2str(plt)],options_);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for jj = 1:nstar0
                    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                end
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(plt));
                fprintf(fidTeX,'\\caption{Smoothed shocks.}');
                fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
            end
        end
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end
    end
    %%
    %%  Smooth observational errors...
    %%
    if options_.prefilter == 1
        yf = atT(bayestopt_.mf,:)+repmat(dataset_info.descriptive.mean',1,gend);
    elseif options_.loglinear == 1
        yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
             trend_coeff*[1:gend];
    else
        yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
             trend_coeff*[1:gend];
    end
    if nvn
        number_of_plots_to_draw = 0;
        index = [];
        for i=1:n_varobs
            if max(abs(measurement_error)) > 0.000000001
                number_of_plots_to_draw = number_of_plots_to_draw + 1;
                index = cat(1,index,i);
            end
            eval(['oo_.SmoothedMeasurementErrors.' options_.varobs{i} ' = measurement_error(i,:)'';']);
        end
        if ~options_.nograph
            [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([M_.fname '_SmoothedObservationErrors.TeX'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation_1.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
                fprintf(fidTeX,' \n');
            end
            for plt = 1:nbplt
                fh = dyn_figure(options_,'Name','Smoothed observation errors');
                NAMES = [];
                if options_.TeX, TeXNAMES = []; end
                nstar0=min(nstar,number_of_plots_to_draw-(nbplt-1)*nstar);
                if gend==1
                    marker_string{1,1}='-ro';
                    marker_string{2,1}='-ko';
                else
                    marker_string{1,1}='-r';
                    marker_string{2,1}='-k';
                end
                for i=1:nstar0
                    k = (plt-1)*nstar+i;
                    subplot(nr,nc,i);
                    plot([1 gend],[0 0],marker_string{1,1},'linewidth',.5)
                    hold on
                    plot(1:gend,measurement_error(index(k),:),marker_string{2,1},'linewidth',1)
                    hold off
                    name = options_.varobs{index(k)};
                    if gend>1
                        xlim([1 gend])
                    end
                    if isempty(NAMES)
                        NAMES = name;
                    else
                        NAMES = char(NAMES,name);
                    end
                    if ~isempty(options_.XTick)
                        set(gca,'XTick',options_.XTick)
                        set(gca,'XTickLabel',options_.XTickLabel)
                    end
                    if options_.TeX
                        idx = strmatch(options_.varobs{index(k)},M_.endo_names,'exact');
                        texname = M_.endo_names_tex(idx,:);
                        if isempty(TeXNAMES)
                            TeXNAMES = ['$ ' deblank(texname) ' $'];
                        else
                            TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                        end
                    end
                    title(name,'Interpreter','none')
                end
                dyn_saveas(fh,[M_.fname '_SmoothedObservationErrors' int2str(plt)],options_);
                if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                    fprintf(fidTeX,'\\begin{figure}[H]\n');
                    for jj = 1:nstar0
                        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                    end
                    fprintf(fidTeX,'\\centering \n');
                    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(plt));
                    fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
                    fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
                    fprintf(fidTeX,'\\end{figure}\n');
                    fprintf(fidTeX,'\n');
                end
            end
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fprintf(fidTeX,'\n');
                fprintf(fidTeX,'%% End of TeX file.\n');
                fclose(fidTeX);
            end
        end
    end
    %%
    %%  Historical and smoothed variabes
    %%
    if ~options_.nograph
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX = fopen([M_.fname '_HistoricalAndSmoothedVariables.TeX'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation_1.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
        fprintf(fidTeX,' \n');
    end
    for plt = 1:nbplt,
        fh = dyn_figure(options_,'Name','Historical and smoothed variables');
        NAMES = [];
        if options_.TeX, TeXNAMES = []; end
        nstar0=min(nstar,n_varobs-(plt-1)*nstar);
        if gend==1
           marker_string{1,1}='-ro';
           marker_string{2,1}='--ko';
        else
           marker_string{1,1}='-r';
           marker_string{2,1}='--k';
        end
        for i=1:nstar0,
            k = (plt-1)*nstar+i;
            subplot(nr,nc,i);
            plot(1:gend,yf(k,:),marker_string{1,1},'linewidth',1)
            hold on
            plot(1:gend,rawdata(:,k),marker_string{2,1},'linewidth',1)
            hold off
            name = options_.varobs{k};
            if isempty(NAMES)
                NAMES = name;
            else
                NAMES = char(NAMES,name);
            end
            if ~isempty(options_.XTick)
                set(gca,'XTick',options_.XTick)
                set(gca,'XTickLabel',options_.XTickLabel)
            end
            if gend>1
                xlim([1 gend])
            end
            if options_.TeX
                idx = strmatch(options_.varobs{k},M_.endo_names,'exact');
                texname = M_.endo_names_tex(idx,:);
                if isempty(TeXNAMES)
                    TeXNAMES = ['$ ' deblank(texname) ' $'];
                else
                    TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                end
            end
            title(name,'Interpreter','none')
        end
        dyn_saveas(fh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)],options_);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:nstar0,
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
            end
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(plt));
            fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
            fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,'\n');
        end
    end
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    end
end

if options_.forecast > 0 && options_.mh_replic == 0 && ~options_.load_mh_file
    dyn_forecast(var_list_,'smoother');
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.fname '_pindx.mat'] ,'pindx');
end

