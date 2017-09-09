function [info,description] = check_posterior_analysis_data(type,M_)
% function [info,description] = check_posterior_analysis_data(type,M_)
% Checks the status of posterior analysis and in particular if files need to be
% created or updated; called by posterior_analysis.m
%
% Inputs:
%   type        [string]        name of the posterior moment considered
%   M_          [structure]     Dynare model structure
%
% Outputs:
%   info        [scalar]        return code
%                                   info = 1; % select_posterior_draws has to be called first.
%                                   info = 2; % _posterior_draws files have to be updated.
%                                   info = 3; % Ok! posterior draws files are up to date ;
%                                   info = 4; % posterior draws have to be processed.
%                                   info = 5; % posterior data files have to be updated.
%                                   info = 6; % Ok (nothing to do ;-)
%   description [string]        Message corresponding to info

% Copyright (C) 2008-2017 Dynare Team
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

info = 0;
if nargout>1
    description = '';
end

[MetropolisFolder, info] = CheckPath('metropolis',M_.dname);

% Get informations about mcmc files.
if info
    disp('check_posterior_analysis_data:: Can''t find any mcmc file!')
    return
end

mhname = get_name_of_the_last_mh_file(M_);
mhdate = get_date_of_a_file([MetropolisFolder filesep mhname]);

% Get informations about _posterior_draws files.
drawsinfo = dir([ MetropolisFolder filesep M_.fname '_posterior_draws*.mat']);
if isempty(drawsinfo)
    info = 1; % select_posterior_draws has to be called first.
    if nargout>1
        description = 'select_posterior_draws has to be called.';
    end
    return
else
    number_of_last_posterior_draws_file = length(drawsinfo);
    pddate = get_date_of_a_file([ MetropolisFolder filesep M_.fname '_posterior_draws' int2str(number_of_last_posterior_draws_file) '.mat']);
    if pddate<mhdate
        info = 2; % _posterior_draws files have to be updated.
        if nargout>1
            description = 'posterior draws files have to be updated.';
        end
        return
    else
        info = 3; % Ok!
        if nargout>1
            description = 'posterior draws files are up to date.';
        end
    end
end

% Get informations about posterior data files.
switch type
  case 'variance'
    generic_post_data_file_name = 'Posterior2ndOrderMoments';
  case 'decomposition'
    generic_post_data_file_name = 'PosteriorVarianceDecomposition';
  case 'correlation'
    generic_post_data_file_name = 'PosteriorCorrelations';
  case 'conditional decomposition'
    generic_post_data_file_name = 'PosteriorConditionalVarianceDecomposition';
  otherwise
    disp('This feature is not yet implemented!')
end
pdfinfo = dir([ MetropolisFolder filesep M_.fname '_' generic_post_data_file_name '*']);
if isempty(pdfinfo)
    info = 4; % posterior draws have to be processed.
    if nargout>1
        description = 'posterior draws files have to be processed.';
    end
    return
else
    number_of_the_last_post_data_file = length(pdfinfo);
    name_of_the_last_post_data_file = ...
        [ pwd filesep MetropolisFolder filesep ...
          M_.fname '_' ...
          generic_post_data_file_name ...
          int2str(number_of_the_last_post_data_file) ...
          '.mat' ];
    pdfdate = get_date_of_a_file(name_of_the_last_post_data_file);
    if pdfdate<pddate
        info = 5; % posterior data files have to be updated.
        if nargout>1
            description = 'posterior data files have to be updated.';
        end
    else
        info = 6; % Ok (nothing to do ;-)
        if nargout>1
            description = 'There is nothing to do';
        end
    end
end