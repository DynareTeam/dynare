function [posterior_sampler_options, options_] = check_posterior_sampler_options(posterior_sampler_options, options_)

% function posterior_sampler_options = check_posterior_sampler_options(posterior_sampler_options, options_)
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

posterior_sampler_options.posterior_sampling_method = options_.posterior_sampling_method;
posterior_sampler_options.proposal_distribution = options_.proposal_distribution;
bounds = posterior_sampler_options.bounds;
invhess = posterior_sampler_options.invhess;

% here are all samplers requiring a proposal distribution
if ~strcmp(posterior_sampler_options.posterior_sampling_method,'slice')
    if ~options_.cova_compute
        error('I Cannot start the MCMC because the Hessian of the posterior kernel at the mode was not computed.')
    end
    if options_.load_mh_file && options_.use_mh_covariance_matrix,
        [junk, invhess] = compute_mh_covariance_matrix;
        posterior_sampler_options.invhess = invhess;
    end
    posterior_sampler_options.parallel_bar_refresh_rate=50;
    posterior_sampler_options.serial_bar_refresh_rate=3;
    posterior_sampler_options.parallel_bar_title='MH';
    posterior_sampler_options.serial_bar_title='Metropolis-Hastings';
    posterior_sampler_options.save_tmp_file=1;
end


% check specific options for slice sampler
if strcmp(options_.posterior_sampling_method,'slice')
    % by default, slice sampler should trigger 
    % mode_compute=0 and
    % mh_replic=100 (much smaller than the default mh_replic=20000 of RWMH)
    % moreover slice must be associated to: 
    %     options_.mh_posterior_mode_estimation = 0;
    % this is done below, but perhaps preprocessing should do this?

    posterior_sampler_options.parallel_bar_refresh_rate=1;
    posterior_sampler_options.serial_bar_refresh_rate=1;
    posterior_sampler_options.parallel_bar_title='SLICE';
    posterior_sampler_options.serial_bar_title='SLICE';
    posterior_sampler_options.save_tmp_file=1;
    posterior_sampler_options = set_default_option(posterior_sampler_options,'rotated',0);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'slice_initialize_with_mode',0);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'use_mh_covariance_matrix',0);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'WR',[]);
    if ~isfield(posterior_sampler_options,'mode'),
        posterior_sampler_options.mode = [];
    else % multimodal case
        posterior_sampler_options.rotated = 1;
    end
    posterior_sampler_options = set_default_option(posterior_sampler_options,'mode_files',[]);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'W1',0.8*(bounds.ub-bounds.lb));

    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
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
                % for multimodal posteriors provide a list of mode files,
                % one per mode. With this info, the code will automatically
                % set the <mode> option. The mode files need only to
                % contain the xparam1 variable.
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
                posterior_sampler_options.W1 = options_list{i,2}*(bounds.ub-bounds.lb);
                
              case 'use_mh_covariance_matrix'
                % in association with <rotated> indicates to use the
                % covariance matrix from previous iterations to define the
                % rotated slice
                % default = 0
                posterior_sampler_options.use_mh_covariance_matrix = options_list{i,2};
                
              otherwise
                warning(['slice_sampler: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    if options_.load_mh_file,
        posterior_sampler_options.slice_initialize_with_mode = 0;
    else
        if ~posterior_sampler_options.slice_initialize_with_mode,
            posterior_sampler_options.invhess=[];
        end
    end    

    if posterior_sampler_options.rotated,
        if isempty(posterior_sampler_options.mode_files) && isempty(posterior_sampler_options.mode), % rotated unimodal
                
            if ~options_.cova_compute && ~(options_.load_mh_file && posterior_sampler_options.use_mh_covariance_matrix)
                skipline()
                disp('I cannot start rotated slice sampler because')
                disp('there is no previous MCMC to load ')
                disp('or the Hessian at the mode is not computed.')
                error('Rotated slice cannot start')
            end
            if isempty(invhess)
                error('oops! This error should not occur, please contact developers.')
            end                
            if options_.load_mh_file && posterior_sampler_options.use_mh_covariance_matrix,
                [junk, invhess] = compute_mh_covariance_matrix;
                posterior_sampler_options.invhess = invhess;
            end
            [V1 D]=eig(invhess);
            posterior_sampler_options.V1=V1;
            posterior_sampler_options.WR=sqrt(diag(D))*3;
        end
    end

    if ~isempty(posterior_sampler_options.mode_files), % multimodal case
        modes = posterior_sampler_options.mode_files; % these can be also mean files from previous parallel slice chains
        for j=1:length( modes ),
            load(modes{j}, 'xparam1')
            mode(j).m=xparam1;
        end
        posterior_sampler_options.mode = mode;
        posterior_sampler_options.rotated = 1;
        posterior_sampler_options.WR=[];
    end
    options_.mh_posterior_mode_estimation = 0;
    
end



