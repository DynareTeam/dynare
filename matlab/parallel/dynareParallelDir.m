function dirlist = dynareParallelDir(filename,PRCDir,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized version of dir() function.
%
% INPUTS
%  o filename   []   ...
%  o PRCDir     []   ...
%  o Parallel   []   ...
%
%  OUTPUTS
%  o dirlist    []   ...
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

dirlist=[];
for indPC=1:length(Parallel),
    if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem), %isunix || (~matlab_ver_less_than('7.4') && ismac),
        if Parallel(indPC).Local==0,
            if ~isempty(Parallel(indPC).Port),
                ssh_token = ['-p ',Parallel(indPC).Port];
            else
                ssh_token = '';
            end
            if exist('OCTAVE_VERSION') % Patch for peculiar behaviour of ssh-ls under Linux.
                % It is necessary to capture the ls warning message.
                % To do it under the ssh protocol it is necessary to redirect the ls message in a text file.
                % The file is 'OctaveStandardOutputMessage.txt' and it is
                % saved in the Model directory.
                [check, ax]=system(['ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' ls ',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',filename, ' 2> OctaveStandardOutputMessage.txt']);
            else
                [check, ax]=system(['ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' ls ',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',filename]);
            end
            if check ~= 0 || ~isempty(strfind(ax,'No such file or directory'));
                ax=[];
            else
                indax=regexp(ax,'\n');
                colax=indax(1);
                rowax=length(indax);
                ax=reshape(ax',[colax rowax])';
                ax=ax(:,1:end-1);
            end
        else

            if exist('OCTAVE_VERSION') % Patch for peculiar behaviour of ls under Linux.
                
                % It is necessary to capture the ls warning message and properly manage the jolly char '*'!
                [check ax]=system(['ls ' ,filename, ' 2> OctaveStandardOutputMessage.txt']);
                
                if check ~= 0 || ~isempty(strfind(ax,'No such file or directory'));
                    ax=[];
                end
            else
                try
                    ax=ls(filename);
                catch
                    ax=[];
                end
            end

        end
    else
        if exist('OCTAVE_VERSION'),     % Patch for peculiar behaviour of ls under Windows.
            if Parallel(indPC).Local==0,
                ax0=dir(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',filename]);
            else
                ax0=dir(filename);
            end
            if isempty(ax0),
                ax='';
            else
                clear ax1;
                for jax=1:length(ax0);
                    ax1{jax}=ax0(jax).name;
                end
                ax=char(ax1{:});
            end

        else
            if Parallel(indPC).Local==0,
                ax=ls(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',filename]);
            else
                ax=ls(filename);
            end
        end
    end
    if isempty(dirlist),
        dirlist=ax;
    elseif ~isempty(ax),
        dirlist = char(dirlist, ax);
    end
end
