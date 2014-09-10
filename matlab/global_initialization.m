function global_initialization()
%function global_initialization()
% initializes global variables and options for DYNARE
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2014 Dynare Team
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

global oo_ M_ options_ estim_params_ bayestopt_ estimation_info ex0_ ys0_

estim_params_ = [];
bayestopt_ = [];
options_.datafile = '';
options_.dataset = [];
options_.verbosity = 1;
options_.terminal_condition = 0;
options_.rplottype = 0;
options_.smpl = 0;
options_.dynatol.f = 1e-5;
options_.dynatol.x = 1e-5;
options_.slowc = 1;
options_.timing = 0;
options_.gstep = ones(2,1);
options_.gstep(1) = 1e-2;
options_.gstep(2) = 1.0;
options_.scalv = 1;
options_.debug = 0;
options_.initval_file = 0;
options_.Schur_vec_tol = 1e-11; % used to find nonstationary variables in Schur decomposition of the
                                % transition matrix
options_.qz_criterium = [];
options_.qz_zero_threshold = 1e-6;
options_.lyapunov_complex_threshold = 1e-15;
options_.solve_tolf = eps^(1/3);
options_.solve_tolx = eps^(2/3);

options_.dp.maxit = 3000;
options_.steady.maxit = 10;
options_.simul.maxit = 10;

options_.mode_check.status = 0;
options_.mode_check.neighbourhood_size = .5;
options_.mode_check.symmetric_plots = 1;
options_.mode_check.number_of_points = 20;
options_.mode_check.nolik = 0;

% Default number of threads for parallelized mex files.
options_.threads.kronecker.A_times_B_kronecker_C = 1;
options_.threads.kronecker.sparse_hessian_times_B_kronecker_C = 1;
options_.threads.local_state_space_iteration_2 = 1;

% steady state
options_.jacobian_flag = 1;

% steady state file
if exist([M_.fname '_steadystate2.m'],'file')
    options_.steadystate_flag = 2;
elseif exist([M_.fname '_steadystate.m'],'file')
    options_.steadystate_flag = 1;
else
    options_.steadystate_flag = 0;
end
options_.steadystate_partial = [];
options_.steadystate.nocheck = 0;

% subset of the estimated deep parameters
options_.ParamSubSet = 'None';

% bvar-dsge
M_.bvar = [];
options_.dsge_var = 0;
options_.dsge_varlag = 4;

% BVAR a la Sims
options_.bvar_replic = 2000;
options_.bvar_prior_tau = 3;
options_.bvar_prior_decay = 0.5;
options_.bvar_prior_lambda = 5;
options_.bvar_prior_mu = 2;
options_.bvar_prior_omega = 1;
options_.bvar_prior_flat = 0;
options_.bvar_prior_train = 0;

% Optimization algorithm [6] gmhmaxlik
gmhmaxlik.iterations = 3;
gmhmaxlik.number = 20000;
gmhmaxlik.nclimb = 200000;
gmhmaxlik.nscale = 200000;
gmhmaxlik.target = 1/3; % Target for the acceptance rate.
options_.gmhmaxlik = gmhmaxlik;

% Request user input.
options_.nointeractive = 0;

% Graphics
options_.graphics.nrows = 3;
options_.graphics.ncols = 3;
options_.graphics.line_types = {'b-'};
options_.graphics.line_width = 1;
options_.graph_format = 'eps';
options_.nodisplay = 0;
options_.nograph = 0;
options_.XTick = [];
options_.XTickLabel = [];
options_.console_mode = 0;
if isoctave
    if sum(get(0,'screensize'))==4
        options_.console_mode = 1;
        options_.nodisplay = 1;
    end
else
    if isunix && (~usejava('jvm') || ~feature('ShowFigureWindows'))
        options_.console_mode = 1;
        options_.nodisplay = 1;
    end
end

% IRFs & other stoch_simul output
options_.irf = 40;
options_.impulse_responses.plot_threshold=1e-10;
options_.relative_irf = 0;
options_.ar = 5;
options_.hp_filter = 0;
options_.hp_ngrid = 512;
options_.nomoments = 0;
options_.nocorr = 0;
options_.periods = 0;
options_.noprint = 0;
options_.SpectralDensity.trigger = 0;
options_.SpectralDensity.plot  = 1; 
options_.SpectralDensity.cutoff  = 150; 
options_.SpectralDensity.sdl = 0.01; 
options_.nofunctions = false;

% Extended path options
%
% Set debug flag
ep.debug = 0;
% Set memory flag
ep.memory = 0;
% Set verbose mode
ep.verbosity = 0;
% Set bytecode flag
ep.use_bytecode = 0;
% Initialization of the perfect foresight equilibrium paths
% * init=0, previous solution is used.
% * init=1, a path generated with the first order reduced form is used.
% * init=2, mix of cases 0 and 1.
ep.init = 0;
% Maximum number of iterations for the deterministic solver.
ep.maxit = 500;
% Number of periods for the perfect foresight model.
ep.periods = 200;
% Default step for increasing the number of periods if needed
ep.step = 50;
% Set check_stability flag
ep.check_stability = 0;
% Define last periods used to test if the solution is stable with respect to an increase in the number of periods.
ep.lp = 5;
% Define first periods used to test if the solution is stable with respect to an increase in the number of periods.
ep.fp = 2;
% Define the distribution for the structural innovations.
ep.innovation_distribution = 'gaussian';
% Set flag for the seed
ep.set_dynare_seed_to_default = 1;
% Set algorithm for the perfect foresight solver
ep.stack_solve_algo = 4;
% Stochastic extended path related options.
ep.stochastic.method = '';
ep.stochastic.algo = 0;
ep.stochastic.quadrature.ortpol = 'hermite';
ep.stochastic.order = 0;
ep.stochastic.quadrature.nodes = 5;
ep.stochastic.quadrature.pruned.status = 0;
ep.stochastic.quadrature.pruned.relative = 1e-5;
ep.stochastic.quadrature.pruned.level = 1e-5;
ep.stochastic.hybrid_order = 0;
% homotopic step in extended path simulations
ep.stochastic.homotopic_steps = true;
% Copy ep structure in options_ global structure
options_.ep = ep;


% Simulations of backward looking models options
%
bnlms.set_dynare_seed_to_default = 1;
bnlms.innovation_distribution = 'gaussian';
options_.bnlms = bnlms;


% Particle filter
%
% Default is that we do not use the non linear kalman filter
particle.status = 0;
% How do we initialize the states?
particle.initialization = 1;
particle.initial_state_prior_std = .1;
% Set the default order of approximation of the model (perturbation).
particle.perturbation = 2;
% Set the default number of particles.
particle.number_of_particles = 5000;
% Set the default approximation order (Smolyak)
particle.smolyak_accuracy = 3;
% By default we don't use pruning
particle.pruning = 0;
% Set the Gaussian approximation method for particles distributions
particle.approximation_method = 'unscented';
% Set unscented parameters alpha, beta and kappa for gaussian approximation
particle.unscented.alpha = .01;
particle.unscented.beta = 2;
particle.unscented.kappa = 0;
% Configuration of resampling in case of particles
particle.resampling.status.systematic = 1;
particle.resampling.status.none = 0;
particle.resampling.status.generic = 0;
particle.resampling.threshold = .5;
particle.resampling.method.kitagawa = 1;
particle.resampling.method.smooth = 0;
particle.resampling.method.stratified = 0;
% Set default algorithm
particle.filter_algorithm.sis = 1;
particle.filter_algorithm.apf = 0;
particle.filter_algorithm.gf = 0;
particle.filter_algorithm.gmf = 0;
particle.filter_algorithm.cpf = 0;
% Approximation of the proposal distribution
particle.proposal_approximation.cubature = 1;
particle.proposal_approximation.unscented = 0;
particle.proposal_approximation.montecarlo = 0;
% Approximation of the particle distribution
particle.distribution_approximation.cubature = 0;
particle.distribution_approximation.unscented = 1;
particle.distribution_approximation.montecarlo = 0;
% Number of partitions for the smoothed resampling method
particle.resampling.number_of_partitions = 200;
% Configuration of the mixture filters
%particle.mixture_method = 'particles' ;
% Size of the different mixtures
particle.mixture_state_variables = 5 ;
particle.mixture_structural_shocks = 1 ;
particle.mixture_measurement_shocks = 1 ;
% Online approach 
particle.liu_west_delta = 0.99 ;
particle.liu_west_chol_sigma_bar = .01 ;
% Copy ep structure in options_ global structure
options_.particle = particle;

% TeX output
options_.TeX = 0;

% Exel
options_.xls_sheet = 1; % Octave does not support the empty string, rather use first sheet
options_.xls_range = '';

% Prior draws
options_.prior_draws = 10000;

options_.forecast = 0;

% Model
options_.linear = 0;

% Deterministic simulation
options_.stack_solve_algo = 0;
options_.markowitz = 0.5;
options_.minimal_solving_periods = 1;
options_.endogenous_terminal_period = 0;
options_.no_homotopy = 0;

% Solution
options_.order = 2;
options_.pruning = 0;
options_.solve_algo = 4;
options_.linear = 0;
options_.replic = 50;
options_.simul_replic = 1;
options_.drop = 100;
% if mjdgges.dll (or .mexw32 or ....) doesn't exist, matlab/qz is added to the path.
% There exists now qz/mjdgges.m that contains the calls to the old Sims code
% Hence, if mjdgges.m is visible exist(...)==2,
% this means that the DLL isn't avaiable and use_qzdiv is set to 1
if exist('mjdgges','file')==2
    options_.use_qzdiv = 1;
else
    options_.use_qzdiv = 0;
end
options_.aim_solver = 0; % i.e. by default do not use G.Anderson's AIM solver, use mjdgges instead
options_.k_order_solver=0; % by default do not use k_order_perturbation but mjdgges
options_.partial_information = 0;
options_.ACES_solver = 0;
options_.conditional_variance_decomposition = [];

% Ramsey policy
options_.ramsey_policy = 0;
options_.timeless = 0;

% estimation
estimation_info.empty_prior = struct(...
    'domain', [], 'interval', [], 'mean', [], ...
    'median', [], 'mode', [], 'shape', [], ...
    'shift', [], 'stdev', [], 'truncate', [], 'variance', []);
estimation_info.empty_options = struct(...
    'bounds',[], 'init', [], 'jscale', []);
estimation_info.subsamples.range = struct('date1', [], 'date2', []);
estimation_info.parameter.prior = estimation_info.empty_prior;
estimation_info.parameter.subsample_prior = estimation_info.empty_prior;
estimation_info.parameter.options = estimation_info.empty_options;
estimation_info.parameter.subsample_options = estimation_info.empty_options;
estimation_info.structural_innovation.prior = estimation_info.empty_prior;
estimation_info.structural_innovation.subsample_prior = estimation_info.empty_prior;
estimation_info.structural_innovation.options = estimation_info.empty_options;
estimation_info.structural_innovation.subsample_options = estimation_info.empty_options;
estimation_info.structural_innovation_corr.prior = estimation_info.empty_prior;
estimation_info.structural_innovation_corr.subsample_prior = estimation_info.empty_prior;
estimation_info.structural_innovation_corr.options = estimation_info.empty_options;
estimation_info.structural_innovation_corr.subsample_options = estimation_info.empty_options;
estimation_info.measurement_error.prior = estimation_info.empty_prior;
estimation_info.measurement_error.subsample_prior = estimation_info.empty_prior;
estimation_info.measurement_error.options = estimation_info.empty_options;
estimation_info.measurement_error.subsample_options = estimation_info.empty_options;
estimation_info.measurement_error_corr.prior = estimation_info.empty_prior;
estimation_info.measurement_error_corr.subsample_prior = estimation_info.empty_prior;
estimation_info.measurement_error_corr.options = estimation_info.empty_options;
estimation_info.measurement_error_corr.subsample_options = estimation_info.empty_options;
estimation_info.subsamples_index = {};
estimation_info.subsamples.range_index = {};
estimation_info.parameter_prior_index = {};
estimation_info.parameter_options_index = {};
estimation_info.parameter.range_index = {};
estimation_info.measurement_error_prior_index = {};
estimation_info.measurement_error_options_index = {};
estimation_info.measurement_error.range_index = {};
estimation_info.structural_innovation_prior_index = {};
estimation_info.structural_innovation_options_index = {};
estimation_info.structural_innovation.range_index = {};
estimation_info.measurement_error_corr_prior_index = {};
estimation_info.measurement_error_corr_options_index = {};
estimation_info.measurement_error_corr.range_index = {};
estimation_info.structural_innovation_corr_prior_index = {};
estimation_info.structural_innovation_corr_options_index = {};
estimation_info.structural_innovation_corr.range_index = {};
options_.initial_period = NaN; %dates(1,1);
options_.dataset.file = [];
options_.dataset.series = [];
options_.dataset.firstobs = dates();
options_.dataset.lastobs = dates();
options_.dataset.nobs = NaN;
options_.dataset.xls_sheet = [];
options_.dataset.xls_range = [];
options_.Harvey_scale_factor = 10;
options_.MaxNumberOfBytes = 1e6;
options_.MaximumNumberOfMegaBytes = 111;
options_.analytic_derivation = 0;
options_.analytic_derivation_mode = 0;
options_.bayesian_irf = 0;
options_.bayesian_th_moments = 0;
options_.diffuse_filter = 0;
options_.filter_step_ahead = [];
options_.filtered_vars = 0;
options_.first_obs = NaN;
options_.nobs = NaN;
options_.kalman_algo = 0;
options_.fast_kalman = 0;
options_.kalman_tol = 1e-10;
options_.use_univariate_filters_if_singularity_is_detected = 1;
options_.riccati_tol = 1e-6;
options_.lik_algo = 1;
options_.lik_init = 1;
options_.load_mh_file = 0;
options_.logdata = 0;
options_.loglinear = 0;
options_.logged_steady_state = 0;
options_.mh_conf_sig = 0.90;
options_.prior_interval = 0.90;
options_.mh_drop = 0.5;
options_.mh_jscale = 0.2;
options_.mh_init_scale = 2*options_.mh_jscale;
options_.mh_mode = 1;
options_.mh_nblck = 2;
options_.mh_recover = 0;
options_.mh_replic = 20000;
options_.recursive_estimation_restart = 0;
options_.MCMC_jumping_covariance='hessian';
options_.use_calibration_initialization = 0;
options_.endo_vars_for_moment_computations_in_estimation=[];

% Prior restrictions
options_.prior_restrictions.status = 0;
options_.prior_restrictions.routine = [];

options_.mode_compute = 4;
options_.mode_file = '';
options_.moments_varendo = 0;
options_.nk = 1;
options_.noconstant = 0;
options_.nodiagnostic = 0;
options_.mh_posterior_mode_estimation = 0;
options_.prefilter = 0;
options_.presample = 0;
options_.prior_trunc = 1e-10;
options_.smoother = 0;
options_.student_degrees_of_freedom = 3;
options_.posterior_max_subsample_draws = 1200;
options_.sub_draws = [];
options_.use_mh_covariance_matrix = 0;
options_.gradient_method = 2;
options_.gradient_epsilon = 1e-6;
options_.posterior_sampling_method = 'random_walk_metropolis_hastings';
options_.proposal_distribution = 'rand_multivariate_normal';
options_.student_degrees_of_freedom = 3;
options_.trace_plot_ma = 200;
options_.mh_autocorrelation_function_size = 30;
options_.plot_priors = 1;
options_.cova_compute = 1;
options_.parallel = 0;
options_.parallel_info.isHybridMatlabOctave = false;
options_.parallel_info.leaveSlaveOpen = 0;
options_.parallel_info.RemoteTmpFolder = '';
options_.number_of_grid_points_for_kde = 2^9;
quarter = 1;
years = [1 2 3 4 5 10 20 30 40 50];
options_.conditional_variance_decomposition_dates = zeros(1,length(years));
for i=1:length(years)
    options_.conditional_variance_decomposition_dates(i) = ...
        (years(i)-1)*4+quarter;
end
options_.filter_covariance = 0;
options_.filter_decomposition = 0;
options_.selected_variables_only = 0;
options_.initialize_estimated_parameters_with_the_prior_mode = 0;
options_.estimation_dll = 0;
% Misc
options_.conf_sig = 0.6;
oo_.exo_simul = [];
oo_.endo_simul = [];
ys0_ = [];
ex0_ = [];
oo_.dr = [];
oo_.exo_steady_state = [];
oo_.exo_det_steady_state = [];
oo_.exo_det_simul = [];

M_.params = [];
M_.endo_histval = [];
M_.Correlation_matrix = [];
M_.Correlation_matrix_ME = [];

% homotopy for steady state
options_.homotopy_mode = 0;
options_.homotopy_steps = 1;
options_.homotopy_force_continue = 0;

% Simplex optimization routine (variation on Nelder Mead algorithm).
simplex.tolerance.x = 1e-4;
simplex.tolerance.f = 1e-4;
simplex.maxiter = 5000;
simplex.maxfcallfactor = 500;
simplex.maxfcall = [];
options_.simplex = simplex;

% CMAES optimization routine.
cmaes.SaveVariables='off';
cmaes.DispFinal='on';
cmaes.WarnOnEqualFunctionValues='no';
cmaes.DispModulo='10';
cmaes.LogModulo='0';
cmaes.LogTime='0';
cmaes.TolFun = 1e-7;
cmaes.TolX = 1e-7;
options_.cmaes = cmaes;

% simpsa optimization routine.
simpsa.TOLFUN = 1e-4;
simpsa.TOLX = 1e-4;
simpsa.TEMP_END = .1;
simpsa.COOL_RATE = 10;
simpsa.INITIAL_ACCEPTANCE_RATIO = .95;
simpsa.MIN_COOLING_FACTOR = .9;
simpsa.MAX_ITER_TEMP_FIRST = 50;
simpsa.MAX_ITER_TEMP_LAST = 2000;
simpsa.MAX_ITER_TEMP = 10;
simpsa.MAX_ITER_TOTAL = 5000;
simpsa.MAX_TIME = 2500;
simpsa.MAX_FUN_EVALS = 20000;
simpsa.DISPLAY = 'iter';
options_.simpsa = simpsa;

% prior analysis
options_.prior_mc = 20000;
options_.prior_analysis_endo_var_list = [];

% did model undergo block decomposition + minimum feedback set computation ?
options_.block = 0;

% model evaluated using a compiled MEX
options_.use_dll = 0;

% model evaluated using bytecode.dll
options_.bytecode = 0;

% if equal to 1 use a fixed point method to solve Sylvester equation (for large scale models)
options_.sylvester_fp = 0;

% convergence criteria to solve iteratively a sylvester equations
options_.sylvester_fixed_point_tol = 1e-12;

% if 1 use a fixed point method to solve Lyapunov equation (for large scale models)
options_.lyapunov_fp = 0;
% if 1 use a doubling algorithm to solve Lyapunov equation (for large scale models)
options_.lyapunov_db = 0;
% if 1 use a square root solver to solve Lyapunov equation (for large scale models)
options_.lyapunov_srs = 0;

% convergence criterion for iteratives methods to solve lyapunov equations
options_.lyapunov_fixed_point_tol = 1e-10;
options_.lyapunov_doubling_tol = 1e-16;

% if equal to 1 use a cycle reduction method to compute the decision rule (for large scale models)
options_.dr_cycle_reduction = 0;

% convergence criterion for iteratives methods to solve the decision rule
options_.dr_cycle_reduction_tol = 1e-7;

% if equal to 1 use a logarithmic reduction method to compute the decision rule (for large scale models)
options_.dr_logarithmic_reduction = 0;

% convergence criterion for iteratives methods to solve the decision rule
options_.dr_logarithmic_reduction_tol = 1e-12;

% convergence criterion for iteratives methods to solve the decision rule
options_.dr_logarithmic_reduction_maxiter = 100;

% dates for historical time series
options_.initial_date.freq = 1;
options_.initial_date.period = 1;
options_.initial_date.subperiod = 0;

% discretionary policy
options_.discretionary_policy = 0;
options_.discretionary_tol = 1e-7;

% Shock decomposition
options_.parameter_set = [];

% Nonlinearfilters
options_.nonlinear_filter = [];

% SBVAR & MS SBVAR initializations:
% SBVAR
options_.ms.vlistlog = [];
options_.ms.restriction_fname = 0;
options_.ms.cross_restrictions = 0;
options_.ms.contemp_reduced_form = 0;
options_.ms.real_pseudo_forecast = 0;
options_.ms.dummy_obs = 0;
options_.ms.ncsk = 0;
options_.ms.indxgforhat = 1;
options_.ms.indxgimfhat = 1;
options_.ms.indxestima = 1;
options_.ms.indxgdls = 1;
options_.ms.cms =0;
options_.ms.ncms = 0;
options_.ms.eq_cms = 1;
options_.ms.banact = 1;
options_.ms.log_var = [];
options_.ms.Qi = [];
options_.ms.Ri = [];
options_.ms.lower_cholesky = 0;
options_.ms.upper_cholesky = 0;
options_.ms.constants_exclusion = 0;
%options_.ms.nstates = 2;
%options_.ms.indxscalesstates = 0;
%options_.ms.q_diag = 0.85;
%options_.ms.flat_prior = 0;
%options_.ms.nstd = 6;
%options_.ms.ninv = 1000;
%options_.ms.indxparr = 1;
%options_.ms.indxovr = 0;
%options_.ms.aband = 1;
%options_.ms.indxap = 1;
%options_.ms.apband = 1;
%options_.ms.indximf = 1;
%options_.ms.imfband = 1;
%options_.ms.indxfore = 0;
%options_.ms.foreband = 0;
%options_.ms.cnum = 0;

% MS SBVAR (and some SBVAR)
options_ = initialize_ms_sbvar_options(M_, options_);

% saved graph formats
options_.graph_save_formats.eps = 1;
options_.graph_save_formats.pdf = 0;
options_.graph_save_formats.fig = 0;

% risky steady state
options_.risky_steadystate = 0;

% endogenous prior
options_.endogenous_prior = 0;
options_.endogenous_prior_restrictions.irf={};
options_.endogenous_prior_restrictions.moment={};

% OSR Optimal Simple Rules
options_.osr.tolf=1e-7;
options_.osr.maxit=1000;

% use GPU
options_.gpu = 0;

%Geweke convergence diagnostics
options_.convergence.geweke.taper_steps=[4 8 15];
options_.convergence.geweke.geweke_interval=[0.2 0.5];

% Options for lmmcp solver
options_.lmmcp = [];

% initialize persistent variables in priordens()
priordens([],[],[],[],[],[],1);
% initialize persistent variables in dyn_first_order_solver()
dyn_first_order_solver();

% Set dynare random generator and seed.
set_dynare_seed('default');


% Create directories
[junk,junk]=mkdir(M_.fname);
[junk,junk]=mkdir([M_.fname '/Output']);

% Load user configuration file.
if isfield(options_, 'global_init_file')
    run(options_.global_init_file);
end

