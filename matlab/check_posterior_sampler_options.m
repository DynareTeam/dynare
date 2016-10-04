function [posterior_sampler_options, options_] = check_posterior_sampler_options(posterior_sampler_options, options_, bounds)

% function [posterior_sampler_options, options_] = check_posterior_sampler_options(posterior_sampler_options, options_, bounds)
% initialization of posterior samplers
%
% INPUTS
%   posterior_sampler_options:       posterior sampler options
%   options_:       structure storing the options

% OUTPUTS
%   posterior_sampler_options:       checked posterior sampler options
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2015 Dynare Team
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


init=0;
if isempty(posterior_sampler_options),
    init=1;
end

if init,
    % set default options and user defined options
    posterior_sampler_options.posterior_sampling_method = options_.posterior_sampler_options.posterior_sampling_method;
    posterior_sampler_options.bounds = bounds;
    
    switch posterior_sampler_options.posterior_sampling_method
        
        case 'random_walk_metropolis_hastings'
            posterior_sampler_options.parallel_bar_refresh_rate=50;
            posterior_sampler_options.serial_bar_refresh_rate=3;
            posterior_sampler_options.parallel_bar_title='RWMH';
            posterior_sampler_options.serial_bar_title='RW Metropolis-Hastings';
            
            % default options
            posterior_sampler_options = add_fields_(posterior_sampler_options,options_.posterior_sampler_options.rwmh);
            
            % user defined options
            if ~isempty(options_.posterior_sampler_options.sampling_opt)
                options_list = read_key_value_string(options_.posterior_sampler_options.sampling_opt);
                for i=1:rows(options_list)
                    switch options_list{i,1}
                        
                        case 'proposal_distribution'
                            if ~(strcmpi(options_list{i,2}, 'rand_multivariate_student') || ...
                                    strcmpi(options_list{i,2}, 'rand_multivariate_normal'))
                                error(['initial_estimation_checks:: the proposal_distribution option to estimation takes either ' ...
                                    'rand_multivariate_student or rand_multivariate_normal as options']);
                            else
                                posterior_sampler_options.proposal_distribution=options_list{i,2};
                            end
                            
                            
                        case 'student_degrees_of_freedom'
                            if options_list{i,2} <= 0
                                error('initial_estimation_checks:: the student_degrees_of_freedom takes a positive integer argument');
                            else
                                posterior_sampler_options.student_degrees_of_freedom=options_list{i,2};
                            end
                            
                        case 'use_mh_covariance_matrix'
                            % indicates to use the covariance matrix from previous iterations to
                            % define the covariance of the proposal distribution
                            % default = 0
                            posterior_sampler_options.use_mh_covariance_matrix = options_list{i,2};
                            options_.use_mh_covariance_matrix = options_list{i,2};
                        case 'scale_file'
                            % load optimal_mh_scale parameter if previous run was with mode_compute=6
                            % will overwrite jscale from set_prior.m
                            if exist(options_list{i,2},'file') || exist([options_list{i,2},'.mat'],'file')
                                tmp = load(options_list{i,2},'Scale');
                                global bayestopt_
                                bayestopt_.mh_jscale = tmp.Scale;
                                options_.mh_jscale = tmp.Scale;
                                bayestopt_.jscale = ones(size(bounds.lb,1),1)*tmp.Scale;
%                                 options_.mh_init_scale = 2*options_.mh_jscale;
                            else
                                error('initial_estimation_checks:: The specified mh_scale_file does not exist.')
                            end
                        case 'save_tmp_file'
                            posterior_sampler_options.save_tmp_file = options_list{i,2};
                        otherwise
                            warning(['rwmh_sampler: Unknown option (' options_list{i,1} ')!'])
                    end
                end
            end
            
        case 'tailored_random_block_metropolis_hastings'
            posterior_sampler_options.parallel_bar_refresh_rate=5;
            posterior_sampler_options.serial_bar_refresh_rate=1;
            posterior_sampler_options.parallel_bar_title='TaRB-MH';
            posterior_sampler_options.serial_bar_title='TaRB Metropolis-Hastings';
            
            % default options
            posterior_sampler_options = add_fields_(posterior_sampler_options,options_.posterior_sampler_options.tarb);
            
            % user defined options
            if ~isempty(options_.posterior_sampler_options.sampling_opt)
                options_list = read_key_value_string(options_.posterior_sampler_options.sampling_opt);
                for i=1:rows(options_list)
                    
                    switch options_list{i,1}
                        
                        case 'proposal_distribution'
                            if ~(strcmpi(options_list{i,2}, 'rand_multivariate_student') || ...
                                    strcmpi(options_list{i,2}, 'rand_multivariate_normal'))
                                error(['initial_estimation_checks:: the proposal_distribution option to estimation takes either ' ...
                                    'rand_multivariate_student or rand_multivariate_normal as options']);
                            else
                                posterior_sampler_options.proposal_distribution=options_list{i,2};
                            end
                            
                            
                        case 'student_degrees_of_freedom'
                            if options_list{i,2} <= 0
                                error('initial_estimation_checks:: the student_degrees_of_freedom takes a positive integer argument');
                            else
                                posterior_sampler_options.student_degrees_of_freedom=options_list{i,2};
                            end
                            
                        case 'mode_compute'
                            posterior_sampler_options.mode_compute=options_list{i,2};
                            
                        case 'optim'
                            posterior_sampler_options.optim_opt=options_list{i,2};

                        case 'new_block_probability'
                            if options_list{i,2}<0 || options_list{i,2}>1
                                error('check_posterior_sampler_options:: The tarb new_block_probability must be between 0 and 1!')
                            else
                                posterior_sampler_options.new_block_probability=options_list{i,2};
                            end
                        case 'scale_file'
                            % load optimal_mh_scale parameter if previous run was with mode_compute=6
                            % will overwrite jscale from set_prior.m
                            if exist(options_list{i,2},'file') || exist([options_list{i,2},'.mat'],'file')
                                tmp = load(options_list{i,2},'Scale');
                                global bayestopt_
                                bayestopt_.mh_jscale = tmp.Scale;
                                options_.mh_jscale = tmp.Scale;
                                bayestopt_.jscale = ones(size(bounds.lb,1),1)*tmp.Scale;
%                                 options_.mh_init_scale = 2*options_.mh_jscale;
                            else
                                error('initial_estimation_checks:: The specified scale_file does not exist.')
                            end
                        case 'save_tmp_file'
                            posterior_sampler_options.save_tmp_file = options_list{i,2};
    
                        otherwise
                            warning(['tarb_sampler: Unknown option (' options_list{i,1} ')!'])
                            
                    end
                    
                end
                
            end
            
        case 'independent_metropolis_hastings'
            posterior_sampler_options.parallel_bar_refresh_rate=50;
            posterior_sampler_options.serial_bar_refresh_rate=3;
            posterior_sampler_options.parallel_bar_title='IMH';
            posterior_sampler_options.serial_bar_title='Ind. Metropolis-Hastings';
            
            % default options
            posterior_sampler_options = add_fields_(posterior_sampler_options,options_.posterior_sampler_options.imh);
            
            % user defined options
            if ~isempty(options_.posterior_sampler_options.sampling_opt)
                options_list = read_key_value_string(options_.posterior_sampler_options.sampling_opt);
                for i=1:rows(options_list)
                    switch options_list{i,1}
                        
                        case 'proposal_distribution'
                            if ~(strcmpi(options_list{i,2}, 'rand_multivariate_student') || ...
                                    strcmpi(options_list{i,2}, 'rand_multivariate_normal'))
                                error(['initial_estimation_checks:: the proposal_distribution option to estimation takes either ' ...
                                    'rand_multivariate_student or rand_multivariate_normal as options']);
                            else
                                posterior_sampler_options.proposal_distribution=options_list{i,2};
                            end
                            
                            
                        case 'student_degrees_of_freedom'
                            if options_list{i,2} <= 0
                                error('initial_estimation_checks:: the student_degrees_of_freedom takes a positive integer argument');
                            else
                                posterior_sampler_options.student_degrees_of_freedom=options_list{i,2};
                            end
                            
                        case 'use_mh_covariance_matrix'
                            % indicates to use the covariance matrix from previous iterations to
                            % define the covariance of the proposal distribution
                            % default = 0
                            posterior_sampler_options.use_mh_covariance_matrix = options_list{i,2};
                            options_.use_mh_covariance_matrix = options_list{i,2};
                            
                        case 'save_tmp_file'
                            posterior_sampler_options.save_tmp_file = options_list{i,2};

                        otherwise
                            warning(['imh_sampler: Unknown option (' options_list{i,1} ')!'])
                    end
                end
            end
            
            
        case 'slice'
            posterior_sampler_options.parallel_bar_refresh_rate=1;
            posterior_sampler_options.serial_bar_refresh_rate=1;
            posterior_sampler_options.parallel_bar_title='SLICE';
            posterior_sampler_options.serial_bar_title='SLICE';
            
            % default options
            posterior_sampler_options = add_fields_(posterior_sampler_options,options_.posterior_sampler_options.slice);
            
            % user defined options
            if ~isempty(options_.posterior_sampler_options.sampling_opt)
                options_list = read_key_value_string(options_.posterior_sampler_options.sampling_opt);
                for i=1:rows(options_list)
                    switch options_list{i,1}
                        case 'rotated'
                            % triggers rotated slice iterations using a covariance
                            % matrix from initial burn-in iterations
                            % must be associated with:
                            % <use_mh_covariance_matrix> or <slice_initialize_with_mode>
                            % default  = 0
                            posterior_sampler_options.rotated = options_list{i,2};
                            
                        case 'mode'
                            % for multimodal posteriors, provide the list of modes as a
                            % matrix, ordered by column, i.e. [x1 x2 x3] for three
                            % modes x1 x2 x3
                            % MR note: not sure this is possible with the
                            % read_key_value_string ???
                            % if this is not possible <mode_files> does to job in any case
                            % This will automatically trigger <rotated>
                            % default = []
                            tmp_mode = options_list{i,2};
                            for j=1:size(tmp_mode,2),
                                posterior_sampler_options.mode(j).m = tmp_mode(:,j);
                            end
                            
                        case 'mode_files'
                            % for multimodal posteriors provide the name of
                            % a file containing a variable array xparams = [nparam * nmodes]
                            % one column per mode. With this info, the code will automatically
                            % set the <mode> option. 
                            % This will automatically trigger <rotated>
                            % default = []
                            posterior_sampler_options.mode_files = options_list{i,2};
                            
                        case 'slice_initialize_with_mode'
                            % the default for slice is to set mode_compute = 0 in the
                            % preprocessor and start the chain(s) from a random location in the prior.
                            % This option first runs the optimizer and then starts the
                            % chain from the mode. Associated with optios <rotated>, it will
                            % use invhess from the mode to perform rotated slice
                            % iterations.
                            % default = 0
                            posterior_sampler_options.slice_initialize_with_mode = options_list{i,2};
                            
                        case 'initial_step_size'
                            % sets the initial size of the interval in the STEPPING-OUT PROCEDURE
                            % the initial_step_size must be a real number in [0, 1],
                            % and it sets the size as a proportion of the prior bounds,
                            % i.e. the size will be initial_step_size*(UB-LB)
                            % slice sampler requires prior_truncation > 0!
                            % default = 0.8
                            if options_list{i,2}<=0 || options_list{i,2}>=1
                                error('check_posterior_sampler_options:: slice initial_step_size must be between 0 and 1')
                            else
                                posterior_sampler_options.initial_step_size=options_list{i,2};
                            end
                        case 'use_mh_covariance_matrix'
                            % in association with <rotated> indicates to use the
                            % covariance matrix from previous iterations to define the
                            % rotated slice
                            % default = 0
                            posterior_sampler_options.use_mh_covariance_matrix = options_list{i,2};
                            options_.use_mh_covariance_matrix = options_list{i,2};

                        case 'save_tmp_file'
                            posterior_sampler_options.save_tmp_file = options_list{i,2};
                            
                        otherwise
                            warning(['slice_sampler: Unknown option (' options_list{i,1} ')!'])
                    end
                end
            end
            
            % slice posterior sampler does not require mode or hessian to run
            % needs to be set to 1 to skip parts in dynare_estimation_1.m
            % requiring posterior maximization/calibrated smoother before MCMC
            options_.mh_posterior_mode_estimation=1;
            
            if ~ posterior_sampler_options.slice_initialize_with_mode
            % by default, slice sampler should trigger
            % mode_compute=0 and
            % mh_replic=100 (much smaller than the default mh_replic=20000 of RWMH)
                options_.mode_compute = 0;
                options_.cova_compute = 0;
            else
                if (isequal(options_.mode_compute,0) && isempty(options_.mode_file) )
                    skipline()
                    disp('check_posterior_sampler_options:: You have specified the option "slice_initialize_with_mode"')
                    disp('check_posterior_sampler_options:: to initialize the slice sampler using mode information')
                    disp('check_posterior_sampler_options:: but no mode file nor posterior maximization is selected,')
                    error('check_posterior_sampler_options:: The option "slice_initialize_with_mode" is inconsistent with mode_compute=0 or empty mode_file.')
                else
                    options_.mh_posterior_mode_estimation=0;
                end
            end
            
            if any(isinf(bounds.lb)) || any(isinf(bounds.ub)),
                skipline()
                disp('some priors are unbounded and prior_trunc is set to zero')
                error('The option "slice" is inconsistent with prior_trunc=0.')
            end
            
            % moreover slice must be associated to:
            %     options_.mh_posterior_mode_estimation = 0;
            % this is done below, but perhaps preprocessing should do this?
            
            if ~isempty(posterior_sampler_options.mode)
                % multimodal case
                posterior_sampler_options.rotated = 1;
                posterior_sampler_options.WR=[];
            end
            %     posterior_sampler_options = set_default_option(posterior_sampler_options,'mode_files',[]);
            
            
            posterior_sampler_options.W1=posterior_sampler_options.initial_step_size*(bounds.ub-bounds.lb);
            if options_.load_mh_file,
                posterior_sampler_options.slice_initialize_with_mode = 0;
            else
                if ~posterior_sampler_options.slice_initialize_with_mode,
                    posterior_sampler_options.invhess=[];
                end
            end
            
            if ~isempty(posterior_sampler_options.mode_files), % multimodal case
                modes = posterior_sampler_options.mode_files; % these can be also mean files from previous parallel slice chains
                load(modes, 'xparams')
                if size(xparams,2)<2,
                    error(['check_posterior_sampler_options:: Variable xparams loaded in file <' modes '> has size [' int2str(size(xparams,1)) 'x' int2str(size(xparams,2)) ']: it must contain at least two columns, to allow multi-modal sampling.'])
                end
                for j=1:size(xparams,2),
                    mode(j).m=xparams(:,j);
                end
                posterior_sampler_options.mode = mode;
                posterior_sampler_options.rotated = 1;
                posterior_sampler_options.WR=[];
            end
            
        otherwise
            error('check_posterior_sampler_options:: Unknown posterior_sampling_method option %s ',posterior_sampler_options.posterior_sampling_method);
    end
    
    return
end

% here are all samplers requiring a proposal distribution
if ~strcmp(posterior_sampler_options.posterior_sampling_method,'slice')
    if ~options_.cova_compute && ~(options_.load_mh_file && posterior_sampler_options.use_mh_covariance_matrix)
        skipline()
        disp('check_posterior_sampler_options:: I cannot start the MCMC because the Hessian of the posterior kernel at the mode was not computed')
        disp('check_posterior_sampler_options:: or there is no previous MCMC to load ')
        error('check_posterior_sampler_options:: MCMC cannot start')
    end
end

if options_.load_mh_file && posterior_sampler_options.use_mh_covariance_matrix,
    [junk, invhess] = compute_mh_covariance_matrix;
    posterior_sampler_options.invhess = invhess;
end



% check specific options for slice sampler
if strcmp(posterior_sampler_options.posterior_sampling_method,'slice')
    invhess = posterior_sampler_options.invhess;
    
    if posterior_sampler_options.rotated,
        if isempty(posterior_sampler_options.mode_files) && isempty(posterior_sampler_options.mode), % rotated unimodal
            
            if ~options_.cova_compute && ~(options_.load_mh_file && posterior_sampler_options.use_mh_covariance_matrix)
                skipline()
                disp('check_posterior_sampler_options:: I cannot start rotated slice sampler because')
                disp('check_posterior_sampler_options:: there is no previous MCMC to load ')
                disp('check_posterior_sampler_options:: or the Hessian at the mode is not computed.')
                error('check_posterior_sampler_options:: Rotated slice cannot start')
            end
            if isempty(invhess)
                error('check_posterior_sampler_options:: This error should not occur, please contact developers.')
            end
            % % %             if options_.load_mh_file && options_.use_mh_covariance_matrix,
            % % %                 [junk, invhess] = compute_mh_covariance_matrix;
            % % %                 posterior_sampler_options.invhess = invhess;
            % % %             end
            [V1, D]=eig(invhess);
            posterior_sampler_options.V1=V1;
            posterior_sampler_options.WR=sqrt(diag(D))*3;
        end
    end
    
    % needs to be re-set to zero otherwise posterior analysis is filtered
    % out in dynare_estimation_1.m
    options_.mh_posterior_mode_estimation = 0;
    
    
else
    
    
end

return

function posterior_sampler_options = add_fields_(posterior_sampler_options, sampler_options)

fnam = fieldnames(sampler_options);
for j=1:length(fnam)
    posterior_sampler_options.(fnam{j}) = sampler_options.(fnam{j});
end



