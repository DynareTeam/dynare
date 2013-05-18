function [results_cell] = execute_posterior_function(functionhandle,M_,options_,oo_,dataset_,estim_params_,bayestopt_,type)
%[results_cell] = execute_posterior_function(functionhandle,M_,options_,oo_,dataset_,estim_params_,bayestopt_,type)% This function executes a given function on draws of the posterior or prior distribution 
%
% INPUTS
%   functionhandle               Handle to the function to be executed
%   M_           [structure]     Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%   options_     [structure]     Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   oo_          [structure]     Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%   dataset_     [structure]     Matlab's structure storing the dataset
%   estim_params_[structure]     Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   bayestopt_   [structure]     Matlab's structure describing the parameter options (initialized by dynare, see @ref{bayestopt_}).
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   results_cell    [cell]     ndrawsx1 cell array storing the results
%                                of the prior/posterior computations

% Copyright (C) 2013 Dynare Team
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

% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
    posterior = 1;
elseif strcmpi(type,'prior')
    DrawsFiles = dir([M_.dname '/prior/draws/' type '_draws*' ]);
    CheckPath('prior/moments',M_.dname);
    posterior = 0;
else
    disp('Execute_posterior_function:: Unknown type!')
    error('');
end
NumberOfDrawsFiles = length(DrawsFiles);

%% initialize output structure
if posterior
    load([M_.dname '/metropolis/' DrawsFiles(1).name ],'pdraws');
else
    load([M_.dname '/prior/draws/' DrawsFiles(1).name ],'pdraws');
end
xparam1=pdraws{1,1};
% get output size
junk=functionhandle(xparam1,M_,options_,oo_,dataset_,estim_params_,bayestopt_);
%initialize cell with number of columns
results_cell=cell(options_.PosteriorSampleSize,size(junk,2));

%% compute function on draws
iter = 1;
for file = 1:NumberOfDrawsFiles
    if posterior
        load([M_.dname '/metropolis/' DrawsFiles(file).name ],'pdraws');
    else
        load([M_.dname '/prior/draws/' DrawsFiles(file).name ],'pdraws');
    end
    NumberOfDraws = rows(pdraws);
    for linee = 1:NumberOfDraws
        M_ = set_all_parameters(pdraws{linee,1},estim_params_,M_);
        [results_cell(iter,:)]=functionhandle(pdraws{linee,1},M_,options_,oo_,dataset_,estim_params_,bayestopt_);
        iter=iter+1;
    end
end

