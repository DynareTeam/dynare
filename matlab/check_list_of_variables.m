function varlist = check_list_of_variables(options_, M_, varlist)
% This function defines, if necessary, the list of endogenous variables
% for which the posterior statistics have to be computed.
%
%
% INPUTS
%
%   options_        [structure]                 Dynare structure.
%   M_              [structure]                 Dynare structure (related to model definition).
%   varlist         [cell of char arrays]       Array of strings with name of the endogenous variables.
%
% OUTPUTS
%   varlist         [cell of char arrays]
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2018 Dynare Team
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

%get uniques

[junk1, junk2, index_uniques] = varlist_indices(varlist, M_.endo_names);
varlist = varlist(index_uniques);

msg = false;
if options_.dsge_var && options_.bayesian_irf
    if ~isempty(varlist)
        for i=1:size(varlist,1)
            idx = strmatch(varlist{i}, options_.varobs, 'exact');
            if isempty(idx)
                disp(sprintf('%s is not an observed variable!', varlist{i}))
                msg = true;
            end
        end
        if ~isequal(size(varlist), length(options_.varobs))
            msg = true;
        end
        if msg
            skipline()
            disp('Posterior IRFs will be computed for all observed variables.')
            skipline()
        end
    end
    varlist = options_.varobs;
    return
end

if ~isempty(varlist) && ~isempty(options_.endo_vars_for_moment_computations_in_estimation)
    error('You cannot use the consider_all_endogenous or consider_all_observed options when listing variables after the estimation command')
elseif isempty(varlist) && ~isempty(options_.endo_vars_for_moment_computations_in_estimation)
    if strcmp(options_.endo_vars_for_moment_computations_in_estimation,'all_endogenous_variables')
        varlist = M_.endo_names(1:M_.orig_endo_nbr);
    elseif strcmp(options_.endo_vars_for_moment_computations_in_estimation,'only_observed_variables')
        varlist = options_.varobs;
    else
        error('Unknown option')
    end
elseif isempty(varlist) && isempty(options_.endo_vars_for_moment_computations_in_estimation)
    skipline()
    disp(['You did not declare endogenous variables after the estimation/calib_smoother command.'])
    cas = '';
    if options_.bayesian_irf
        cas = 'Posterior IRFs';
    end
    if options_.moments_varendo
        if isempty(cas)
            cas = 'Posterior moments';
        else
            cas = [cas, ', posterior moments'];
        end
    end
    if options_.smoother
        if isempty(cas)
            cas = 'Smoothed variables';
        else
            cas = [cas, ', smoothed variables'];
        end
    end
    if ~isempty(options_.filter_step_ahead)
        if isempty(cas)
            cas = 'k-step ahead filtered variables';
        else
            cas = [cas, ', k-step ahead filtered variables'];
        end
    end
    if options_.forecast
        if isempty(cas)
            cas = 'Forecasts';
        else
            cas = [cas, ' and forecasts'];
        end
    end
    if ~isempty(cas)
        str = sprintf('%s will be computed for the %s endogenous variables of your model', cas, num2str(M_.orig_endo_nbr));
        str = sprintf('%s, this can take a long time ....', str);
        format_text(str, 10)
        if options_.nointeractive
            % Default behaviour is to consider all the endogenous variables.
            varlist = M_.endo_names(1:M_.orig_endo_nbr);
        else
            choice = [];
            while isempty(choice)
                skipline(2)
                disp('Choose one of the following options:')
                skipline()
                disp(' [1] Consider all the endogenous variables.')
                disp(' [2] Consider all the observed endogenous variables.')
                disp(' [3] Stop Dynare and change the mod file.')
                skipline()
                choice = input('options [default is 1] =  ');
                if isempty(choice)
                    choice=1;
                end
                if choice==1
                    varlist = M_.endo_names(1:M_.orig_endo_nbr);
                elseif choice==2
                    varlist = options_.varobs;
                elseif choice==3
                    varlist = NaN;
                else
                    skipline()
                    disp('YOU HAVE TO ANSWER 1, 2 or 3!')
                    skipline()
                end
            end
        end
        if ~iscell(varlist) && isnan(varlist)
            edit([M_.fname '.mod'])
        end
        skipline()
    end
end


function format_text(remain, max_number_of_words_per_line)
index = 0;
line_of_text = [];
while ~isempty(remain)
    [token, remain] = strtok(remain);
    index = index+1;
    if isempty(line_of_text)
        line_of_text = token;
    else
        line_of_text = [line_of_text , ' ' , token];
    end
    if index==max_number_of_words_per_line
        disp(line_of_text)
        index = 0;
        line_of_text = [];
    end
end
if index<max_number_of_words_per_line
    disp(line_of_text)
end