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
        invhess = compute_mh_covariance_matrix;
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
    posterior_sampler_options.parallel_bar_refresh_rate=1;
    posterior_sampler_options.serial_bar_refresh_rate=1;
    posterior_sampler_options.parallel_bar_title='SLICE';
    posterior_sampler_options.serial_bar_title='SLICE';
    posterior_sampler_options.save_tmp_file=1;
    posterior_sampler_options = set_default_option(posterior_sampler_options,'rotated',0);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'slice_initialize_with_mode',0);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'use_slice_covariance_matrix',0);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'WR',[]);
    if ~isfield(posterior_sampler_options,'mode'),
        posterior_sampler_options.mode = [];
    else % multimodal case
        posterior_sampler_options.rotated = 1;
    end
    posterior_sampler_options = set_default_option(posterior_sampler_options,'mode_files',[]);
    posterior_sampler_options = set_default_option(posterior_sampler_options,'W1',0.8*(bounds.ub-bounds.lb));
   
    if options_.load_mh_file,
        posterior_sampler_options.slice_initialize_with_mode = 0;
    end    

    if posterior_sampler_options.rotated,
        if isempty(posterior_sampler_options.mode_files) && ~isfield(posterior_sampler_options,'mode'), % rotated unimodal
                
            if ~options_.cova_compute && ~(options_.load_mh_file && options_.use_mh_covariance_matrix)
                skipline()
                disp('I cannot start rotated slice sampler because')
                disp('there is no previous MCMC to load ')
                disp('or the Hessian at the mode is not computed.')
                error('Rotated slice cannot start')
            end
            if isempty(invhess)
                error('oops! This error should not occur, please contact developers.')
            end                
            if options_.load_mh_file && posterior_sampler_options.use_slice_covariance_matrix,
                invhess = compute_mh_covariance_matrix;
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



