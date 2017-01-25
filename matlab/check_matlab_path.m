function check_matlab_path(change_path_flag)
    
% Copyright (C) 2015-2016 Dynare Team
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

if ~nargin || isempty(change_path_flag)
    change_path_flag = true;
end

% Get path to dynare/matlab folder.
DYNARE_PATH = strrep(which('dynare'),'dynare.m','');

if isempty(DYNARE_PATH)
    % Nothing to do here (this case should not happen)
    disp('dynare.m is not in the Matlab''s path.')
    return
else
    % Removes trailing slash.
    DYNARE_PATH = DYNARE_PATH(1:end-1);
end

% Get matlab path
MATLAB_PATH = path();
if isoctave
    % Octave always has '.' at the top of the path, so remove '.:'
    MATLAB_PATH = MATLAB_PATH(3:end);
end

% Position of DYNARE_PATH in MATLAB_PATH
idDYNARE = strfind(MATLAB_PATH,DYNARE_PATH);

if isempty(idDYNARE)
    disp('dynare.m is not in the Matlab''s path.')
    return
else
    if isequal(length(idDYNARE),1)
        if isequal(idDYNARE, 1)
            % Dynare is on top of matlab's path! Nothing to do here...
            return
        else
            str0 = sprintf('Dynare is not on top of matlab''s path!');
            % Check that this will not create a problem
            MATLAB_PATH_ = path2cell(MATLAB_PATH);
            DYNARE_ROUTINES = getallroutinenames(DYNARE_PATH, getalldirectories(DYNARE_PATH));
            MATLAB_ROUTINES = {};
            for i=1:position(idDYNARE, MATLAB_PATH)
                TMP_MATLAB_ROUTINES = getallroutinenames(MATLAB_PATH_{i});
                MATLAB_ROUTINES = { MATLAB_ROUTINES{:}  TMP_MATLAB_ROUTINES{:} };
            end
            COMMON_ROUTINES = intersect(MATLAB_ROUTINES, DYNARE_ROUTINES);
            if ~isempty(COMMON_ROUTINES)
                warning off backtrace
                skipline()
                if length(COMMON_ROUTINES)==1
                    warning(sprintf('%s This can cause problems because the Dynare version of %s will be overriden.', str0, COMMON_ROUTINES{1}));
                else
                    str1 = repmat('%s, ', 1, length(COMMON_ROUTINES)-1);
                    str2 = 'and %s ';
                    str3 = sprintf(['%s This can cause problems because the Dynare versions of ' str1, str2, 'will be overriden.'], str0, COMMON_ROUTINES{:});
                    warning(str3);
                end
                if change_path_flag
                    skipline()
                    msg = sprintf('I put %s on top of your matlab''s path. Note that this is a', DYNARE_PATH);
                    msg = sprintf(' %s a temporary change (ie will not affect future matlab''s session).', msg);
                    msg = sprintf(' %s If the ordering was intentional, ie if you really want to override the routines distributed with Dynare,', msg);
                    msg = sprintf(' %s you can change this behaviour using option nopathchange (see the reference manual).', msg);
                    warning(msg);
                    skipline()
                    rmpath(DYNARE_PATH)
                    addpath(DYNARE_PATH)    
                end
                warning on backtrace
            end
        end
    else
        % Check that the user did not put all the subfolders in the path.
        % => If DYNARE_PATH/qz is in the path while mjdgges dll is available
        % it most likely means that user wrongly put all subfolders in the 
        % matlab's path!
        mexpath = add_path_to_mex_files([DYNARE_PATH filesep], false);
        MATLAB_PATH = path2cell(MATLAB_PATH);
        for i=1:length(mexpath)
            if exist([mexpath{i} filesep 'mjdgges.' mexext],'file') && ismember([DYNARE_PATH filesep 'qz'],MATLAB_PATH)
                msg = sprintf(['You  put all the dynare/matlab subfolders in matlab''s path! Only ' ...
                               'the dynare/matlab folder (without subfolders)\nshould be in the ' ...
                               'path, Dynare will automatically add any required subfolders in the ' ...
                               'path.']);
                error(msg)
            end
        end
    end
end

function q = path2cell(p)
    % Converts the output of path() to a cell 
    s = strfind(p,pathsep);
    n = length(s)+1;
    q = cell(n,1);
    q(1) = {p(1:s(1)-1)};
    q(n) = {p(s(end)+1:end)};
    for i=2:n-1
        q(i) = {p(s(i-1)+1:s(i)-1)}; 
    end
    
function flist = getallroutinenames(p, excludedsubfolders)
    if nargin<2
        excludedsubfolders = {};
    end
    flist={};
    %get m-files in this directory
    dd = dir([p,filesep '*.m']);
    temp=struct2cell(dd);
    flist=[flist temp(1,:)];
    %deal with subdirectories
    dlist=getalldirectories(p,excludedsubfolders); %first call with excluded directories
    for ii=1:length(dlist)
        flist=[flist getallroutinenames([ p filesep dlist{ii}])]; %recursive calls without subfolders
    end
    

function dlist = getalldirectories(p,excluded_directories)
    if nargin<2
        excluded_directories = {};
    end
    dd = dir(p);
    dir_result=struct2cell(dd);
    directory_indicator=cell2mat(dir_result(4,:));
    dlist = dir_result(1,directory_indicator==1 & ~strcmp('.',dir_result(1,:)) & ~strcmp('..',dir_result(1,:)) & ~ismember(dir_result(1,:),excluded_directories));

function n = position(i, currentpath)
    n = length(strfind(currentpath(1:i), pathsep));