function oo_=execute_prior_posterior_function(posterior_function_name,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info,type)
%[oo_] = execute_prior_posterior_function(posterior_function_name,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info,type)
% This function executes a given function on draws of the posterior or prior distribution
%
% INPUTS
%   functionhandle               Handle to the function to be executed
%   M_           [structure]     Matlab/Octave structure describing the Model (initialized by dynare, see @ref{M_}).
%   options_     [structure]     Matlab/Octave structure describing the options (initialized by dynare, see @ref{options_}).
%   oo_          [structure]     Matlab/Octave structure gathering the results (initialized by dynare, see @ref{oo_}).
%   estim_params_[structure]     Matlab/Octave structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   bayestopt_   [structure]     Matlab/Octave structure describing the parameter options (initialized by dynare, see @ref{bayestopt_}).
%   dataset_     [structure]     Matlab/Octave structure storing the dataset
%   dataset_info [structure]     Matlab/Octave structure storing the information about the dataset
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   oo_          [structure]     Matlab/Octave structure gathering the results (initialized by dynare, see @ref{oo_}).

% Copyright (C) 2013-2015 Dynare Team
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

[directory,basename,extension] = fileparts(posterior_function_name);
if isempty(extension)
    extension = '.m';
end
fullname = [basename extension];
if ~strcmp(extension,'.m') %if not m-file
    error('The Posterior Function is not an m-file.')
elseif ~exist(fullname,'file') %if m-file, but does not exist
    error(['The Posterior Function ', fullname ,' was not found. Check the spelling.']);
end
%Create function handle
functionhandle=str2func(posterior_function_name);

prior = true;
n_draws=options_.sampling_draws;
% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    %% discard first mh_drop percent of the draws:
    CutSample(M_, options_, estim_params_);
    %% initialize metropolis draws
    options_.sub_draws=n_draws; %set draws for sampling; changed value is not returned to base workspace
    [error_flag,junk,options_]= metropolis_draw(1,options_,estim_params_,M_);
    if error_flag
        error('EXECUTE_POSTERIOR_FUNCTION: The draws could not be initialized')
    end
    n_draws=options_.sub_draws;
    prior = false;
elseif strcmpi(type,'prior')
    prior_draw(bayestopt_, options_.prior_trunc);
else
    error('EXECUTE_POSTERIOR_FUNCTION: Unknown type!')
end

%get draws for later use
first_draw=GetOneDraw(type,M_,estim_params_,oo_,options_,bayestopt_);
parameter_mat=NaN(n_draws,length(first_draw));
parameter_mat(1,:)=first_draw;
for draw_iter=2:n_draws
    parameter_mat(draw_iter,:) = GetOneDraw(type,M_,estim_params_,oo_,options_,bayestopt_);
end

% get output size
try
    junk=functionhandle(parameter_mat(1,:),M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info);
catch err
    fprintf('\nEXECUTE_POSTERIOR_FUNCTION: Execution of prior/posterior function led to an error. Execution cancelled.\n')
    rethrow(err)
end

%initialize cell with number of columns
results_cell=cell(n_draws,size(junk,2));

%% compute function on draws
for draw_iter = 1:n_draws
    M_ = set_all_parameters(parameter_mat(draw_iter,:),estim_params_,M_);
    [results_cell(draw_iter,:)]=functionhandle(parameter_mat(draw_iter,:),M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info);
end

if prior
    oo_.prior_function_results = results_cell;
else
    oo_.posterior_function_results = results_cell;
end
