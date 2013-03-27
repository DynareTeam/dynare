function fileList = dynareParallelListAllFiles(dirName,PRCDir,Parallel)

% PARALLEL CONTEXT
% In a parallel context, this function searches recursively through all subdirectories
% of the given directory 'dirName' and then return a list of all file
% finds in 'dirName'.
%
%
% INPUTS
%  o dirName           []   ...
%  o PRCDir            []   ...
%  o Parallel          []   ...
%
%  OUTPUTS
%  o fileList          []   ...
%
% Copyright (C) 2009-2011 Dynare Team
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



if (~ispc || strcmpi('unix',Parallel.OperatingSystem))
   
    fileList={};
    currentPath=[];
    if ~isempty(Parallel.Port),
        ssh_token = ['-p ',Parallel.Port];
    else
        ssh_token = '';
    end

     % Get the data for the current remote directory.
    [Flag fL]=system(['ssh ',ssh_token,' ',' ',Parallel.UserName,'@',Parallel.ComputerName,' ls ',Parallel.RemoteDirectory,'/',PRCDir, ' -R -p -1']);

    % Format the return value fL.
    
    nL=regexp(fL,'\n');
    start=1;
    j=1;

    for (i=1:length(nL))
        
        stringProcessed=fL(start:nL(i)-1);
        
        if isempty(stringProcessed)
            start=nL(i)+1;
            continue
        end
            
        if strfind(stringProcessed,'/')
            if strfind(stringProcessed,PRCDir)
                DD=strfind(stringProcessed,':');
                stringProcessed(DD)='/';
                currentPath=stringProcessed;
            end
            start=nL(i)+1;
            continue
        end
        
        fileList{j,1}=[currentPath stringProcessed];
        start=nL(i)+1;
        j=j+1;
    end 


else
    if (strmatch(dirName, 'Root')==1) % First call in Windows enviroment
        dirName=(['\\',Parallel.ComputerName,'\',Parallel.RemoteDrive,'$\',Parallel.RemoteDirectory,'\',PRCDir]);
    end
    % Get the data for the current directory and exstract from it the index
    % for directories:
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];

    % Get a list of the files:
    fileList = {dirData(~dirIndex).name}';

    % Build the path files:
    if ~isempty(fileList)
        fileList = cellfun(@(x) fullfile(dirName,x),...
            fileList,'UniformOutput',false);
    end

    % Get a list of the subdirectories:
    subDirs = {dirData(dirIndex).name};

    % Find index of subdirectories that are not '.' or '..':
    validIndex = ~ismember(subDirs,{'.','..'});

    % Loop over valid subdirectories (i.e. all subdirectory without '.' and
    % '..':

    for iDir = find(validIndex)
        % Get the subdirectory path:
        nextDir = fullfile(dirName,subDirs{iDir});

        % Recursively call dynareParallelListAllFiles:
        fileList = [fileList; dynareParallelListAllFiles(nextDir,PRCDir,Parallel)];
    end
end


