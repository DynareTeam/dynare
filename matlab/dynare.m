function dynare(fname, varargin)
%       This command runs dynare with specified model file in argument
%       Filename.
%       The name of model file begins with an alphabetic character, 
%       and has a filename extension of .mod or .dyn.
%       When extension is omitted, a model file with .mod extension
%       is processed.
%
% INPUTS
%   fname:      file name
%   varargin:   list of arguments following fname
%             
% OUTPUTS
%   none
%        
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2013 Dynare Team
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

if strcmpi(fname,'help')
    skipline()
    disp(['This is dynare version ' dynare_version() '.'])
    skipline()
    disp('USAGE: dynare FILENAME[.mod,.dyn] [OPTIONS]')
    skipline()
    disp('dynare executes instruction included in FILENAME.mod.')
    disp('See the reference manual for the available options.')
    return
end

warning_config()

if exist('OCTAVE_VERSION')
    if octave_ver_less_than('3.6.0')
        warning('This version of Dynare has only been tested on Octave 3.6.0 and above. Since your Octave version is older than that, Dynare may fail to run, or give unexpected results. Consider upgrading your Octave installation.');
    end
else
    if matlab_ver_less_than('7.3')
        warning('This version of Dynare has only been tested on MATLAB 7.3 (R2006b) and above. Since your MATLAB version is older than that, Dynare may fail to run, or give unexpected results. Consider upgrading your MATLAB installation, or switch to Octave.');
    end
end

% disable output paging (it is on by default on Octave)
more off

% sets default format for save() command
if exist('OCTAVE_VERSION')
    default_save_options('-mat')
end

% detect if MEX files are present; if not, use alternative M-files
dynareroot = dynare_config;

if nargin < 1
    error('DYNARE: you must provide the name of the MOD file in argument')
end

if ~ischar(fname)
    error('DYNARE: argument of dynare must be a text string')
end

% Testing if file have extension
% If no extension default .mod is added
if isempty(strfind(fname,'.'))
    fname1 = [fname '.dyn'];
    d = dir(fname1);
    if length(d) == 0
        fname1 = [fname '.mod'];
    end
    fname = fname1;
    % Checking file extension
else
    if ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.MOD') ...
            && ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.DYN')
        error('DYNARE: argument must be a filename with .mod or .dyn extension')
    end;
end;
d = dir(fname);
if length(d) == 0
    fprintf('\nThe file %s could not be located in the "Current Folder". Check whether you typed in the correct filename\n',fname)
    fprintf('and whether the file is really located in the "Current Folder".\n')
    error(['DYNARE: can''t open ' fname])
elseif ~isempty(strfind(fname,'\')) || ~isempty(strfind(fname,'/'))
    fprintf('\nIt seems you are trying to call a mod-file not located in the "Current Folder". This is not possible.\n')
    fprintf('Please set your "Current Folder" to the folder where the mod-file is located.\n')
    error(['DYNARE: can''t open ' fname, '. It seems to be located in a different folder (or has an invalid filename).'])        
end

% pre-dynare-preprocessor-hook
if exist([fname(1:end-4) '_pre_dynare_preprocessor_hook.m'],'file')
    eval([fname(1:end-4) '_pre_dynare_preprocessor_hook'])
end

command = ['"' dynareroot 'dynare_m" ' fname] ;
for i=2:nargin
    command = [command ' ' varargin{i-1}];
end

[status, result] = system(command);
disp(result)
if ismember('onlymacro', varargin)
    disp('Preprocesser stopped after macroprocessing step because of ''onlymacro'' option.');
    return;
end

% post-dynare-prerocessor-hook
if exist([fname(1:end-4) '_post_dynare_preprocessor_hook.m'],'file')
    eval([fname(1:end-4) '_post_dynare_preprocessor_hook'])
end

% Save preprocessor result in logfile (if `no_log' option not present)
no_log = 0;
for i=2:nargin
    no_log = no_log || strcmp(varargin{i-1}, 'nolog');
end
if ~no_log
    logname = [fname(1:end-4) '.log'];
    fid = fopen(logname, 'w');
    fprintf(fid, '%s', result);
    fclose(fid);
end

if status
    % Should not use "error(result)" since message will be truncated if too long
    error('DYNARE: preprocessing failed')
end

if ~ isempty(find(abs(fname) == 46))
    fname = fname(:,1:find(abs(fname) == 46)-1) ;
end
evalin('base',fname) ;
