function dynareroot = dynare_config(path_to_dynare,verbose)
%function dynareroot = dynare_config(path_to_dynare)
% This function tests the existence of valid mex files (for qz
% decomposition, solution to sylvester equation and kronecker
% products...) and, if needed, add paths to the matlab versions
% of these routines.
% Also adds other directories to the path.
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2014 Dynare Team
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

if nargin && ~isempty(path_to_dynare)
    addpath(path_to_dynare);
end
dynareroot = strrep(which('dynare'),'dynare.m','');

origin = pwd();
cd([dynareroot '/..'])

if ~nargin || nargin==1
    verbose = 1;
end


addpath([dynareroot '/distributions/'])
addpath([dynareroot '/kalman/'])
addpath([dynareroot '/kalman/likelihood'])
addpath([dynareroot '/AIM/'])
addpath([dynareroot '/partial_information/'])
addpath([dynareroot '/ms-sbvar/'])
addpath([dynareroot '/ms-sbvar/identification/'])
addpath([dynareroot '../contrib/ms-sbvar/TZcode/MatlabFiles/'])
addpath([dynareroot '/parallel/'])
addpath([dynareroot '/particle/'])
addpath([dynareroot '/gsa/'])
addpath([dynareroot '/ep/'])
addpath([dynareroot '/lmmcp/'])
addpath([dynareroot '/utilities/doc/'])
addpath([dynareroot '/utilities/tests/src/'])
addpath([dynareroot '/utilities/dates/'])
addpath([dynareroot '/utilities/dataset/'])
addpath([dynareroot '/utilities/general/'])
addpath([dynareroot '/utilities/dseries/'])
addpath([dynareroot '/reports/'])

% For functions that exist only under some Octave versions
% or some MATLAB versions, and for which we provide some replacement functions

if ~isoctave
    % Replacements for rows(), columns() and issquare() (inexistent under MATLAB)
    addpath([dynareroot '/missing/rows_columns'])
    addpath([dynareroot '/missing/issquare'])
    % Replacement for vec() (inexistent under MATLAB)
    addpath([dynareroot '/missing/vec'])
    if ~user_has_matlab_license('statistics_toolbox')
        % Replacements for functions of the stats toolbox
        addpath([dynareroot '/missing/stats/'])
    end
end

% ordeig() doesn't exist in Octave
if isoctave
    addpath([dynareroot '/missing/ordeig'])
end

% ilu is missing in Octave
if isoctave
    addpath([dynareroot '/missing/ilu'])
end

% strjoin is missing in older versions of MATLAB and in Octave < 3.8
if (isoctave && octave_ver_less_than('3.8')) || ...
        (~isoctave && matlab_ver_less_than('8.1'))
    addpath([dynareroot '/missing/strjoin'])
end

% nanmean is in Octave Forge Statistics package and in MATLAB Statistics
% toolbox
if (isoctave && ~user_has_octave_forge_package('statistics')) ...
    || (~isoctave && ~user_has_matlab_license('statistics_toolbox'))
    addpath([dynareroot '/missing/nanmean'])
end

% Add path to MEX files
if isoctave
    addpath([dynareroot '../mex/octave/']);
else
    % Add win32 specific paths for Dynare Windows package
    if strcmp(computer, 'PCWIN')
        mexpath = [dynareroot '../mex/matlab/win32-7.5-8.4'];
        if exist(mexpath, 'dir')
            addpath(mexpath)
        end
    end

    % Add win64 specific paths for Dynare Windows package
    if strcmp(computer, 'PCWIN64')
        if matlab_ver_less_than('7.8')
            mexpath = [dynareroot '../mex/matlab/win64-7.5-7.7'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        else
            mexpath = [dynareroot '../mex/matlab/win64-7.8-8.4'];
            if exist(mexpath, 'dir')
                addpath(mexpath)
            end
        end
    end

    if strcmp(computer, 'MACI')
        mexpath = [dynareroot '../mex/matlab/osx32-7.5-7.11'];
        if exist(mexpath, 'dir')
            addpath(mexpath)
        end
    end

    if strcmp(computer, 'MACI64')
        mexpath = [dynareroot '../mex/matlab/osx64'];
        if exist(mexpath, 'dir')
            addpath(mexpath)
        end
    end

    % Add generic MATLAB path (with higher priority than the previous ones)
    addpath([dynareroot '../mex/matlab/']);
end

%% Set mex routine names
mex_status = cell(1,3);
mex_status(1,1) = {'mjdgges'};
mex_status(1,2) = {'qz'};
mex_status(1,3) = {'Generalized QZ'};
mex_status(2,1) = {'gensylv'};
mex_status(2,2) = {'gensylv'};
mex_status(2,3) = {'Sylvester equation solution'};
mex_status(3,1) = {'A_times_B_kronecker_C'};
mex_status(3,2) = {'kronecker'};
mex_status(3,3) = {'Kronecker products'};
mex_status(4,1) = {'sparse_hessian_times_B_kronecker_C'};
mex_status(4,2) = {'kronecker'};
mex_status(4,3) = {'Sparse kronecker products'};
mex_status(5,1) = {'local_state_space_iteration_2'};
mex_status(5,2) = {'particle/local_state_space_iteration'};
mex_status(5,3) = {'Local state space iteration (second order)'};
number_of_mex_files = size(mex_status,1);
%% Remove some directories from matlab's path. This is necessary if the user has
%% added dynare_v4/matlab with the subfolders. Matlab has to ignore these
%% subfolders if valid mex files exist.
matlab_path = path;
for i=1:number_of_mex_files
    test = strfind(matlab_path,[dynareroot mex_status{i,2}]);
    action = length(test);
    if action
        rmpath([dynareroot mex_status{i,2}]);
        matlab_path = path;
    end
end
%% Test if valid mex files are available, if a mex file is not available
%% a matlab version of the routine is included in the path.
if verbose
    skipline()
    disp('Configuring Dynare ...')
end

for i=1:number_of_mex_files
    test = (exist(mex_status{i,1},'file') == 3);
    if ~test
        addpath([dynareroot mex_status{i,2}]);
        message = '[m]   ';
    else
        message = '[mex] ';
    end
    if verbose
        disp([ message mex_status{i,3} '.' ])
    end
end

% Test if bytecode DLL is present
if exist('bytecode', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'Bytecode evaluation.' ])
end

% Test if k-order perturbation DLL is present
if exist('k_order_perturbation', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'k-order perturbation solver.' ])
end

% Test if dynare_simul_ DLL is present
if exist('dynare_simul_', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'k-order solution simulation.' ])
end

% Test if qmc_sequence DLL is present
if exist('qmc_sequence', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'Quasi Monte-Carlo sequence (Sobol).' ])
end

% Test if MS-SBVAR DLL is present
if exist('ms_sbvar_command_line', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'Markov Switching SBVAR.' ])
    skipline()
end

% Save empty dates and dseries objects (necessary if a mod file is not preprocessed).
dates('initialize');
dseries('initialize');

cd(origin);
