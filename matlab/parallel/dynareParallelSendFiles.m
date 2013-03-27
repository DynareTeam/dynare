function dynareParallelSendFiles(NamFileInput,PRCDir,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized mono-directional (Local to Remote) version of copy()
% function.
%
%
% INPUTS
%  o NamFileInput   []   ...
%  o PRCDir         []   ...
%  o Parallel       []   ...
%
%  OUTPUTS
%  None
%
%
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


if ischar(NamFileInput),
    for j=1:size(NamFileInput,1),
        NamFile(j,:)={'',deblank(NamFileInput(j,:))};
    end
    NamFileInput = NamFile;
end

for indPC=1:length(Parallel),
    if Parallel(indPC).Local==0,
        if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem), %isunix || (~matlab_ver_less_than('7.4') && ismac),
            if ~isempty(Parallel(indPC).Port),
                ssh_token = ['-p ',Parallel(indPC).Port];
            else
                ssh_token = '';
            end
            if ~isempty(Parallel(indPC).Port),
                scp_token = ['-P ',Parallel(indPC).Port];
            else
                scp_token = '';
            end
            for jfil=1:size(NamFileInput,1),
                if ~isempty(NamFileInput{jfil,1})
                    [NonServeL NonServeR]=system(['ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' mkdir -p ',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1}]);
                end
                [NonServeL NonServeR]=system(['scp ',scp_token,' ',NamFileInput{jfil,1},NamFileInput{jfil,2},' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,':',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1}]);
            end
        else
            for jfil=1:size(NamFileInput,1)
                if ~isempty(NamFileInput{jfil,1})
                    if isempty(dynareParallelDir(NamFileInput{jfil,1},PRCDir,Parallel(indPC)))
                        
                        if exist('OCTAVE_VERSION') % Patch for peculiar behaviour of mkdir under Windows.
                            
                            % It is Necessary because Octave is not able to
                            % create two nested directory at the same time.
                            
                           % Remove (if present) the '/' chars. Can be easily transformed
                           % in a function.
                           
                            NamFileInputTemp=NamFileInput{jfil,1};
                            while(1)
                                Bs=strfind(NamFileInputTemp,'/');
                                if isempty(Bs),
                                    break;
                                else
                                    NamFileInputTemp(1,Bs)='\';
                                end
                            end
                            
                            [NonServeL NonServeR]=system(['mkdir \\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInputTemp]);

                        else
                            mkdir(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInput{jfil,1}]);
                        end
                    end
                end
                
                if exist('OCTAVE_VERSION') % Patch for peculiar behaviour copyfile ls under Windows.
                    
                    % It is Necessary because Octave is not able to
                    % use the jolly char '*' with copyfile.
                    
                      % Remove (if present) the '/' chars. Can be easily transformed
                      % in a function.
                    
                    NamFileInputTemp=NamFileInput{jfil,1};
                    while(1)
                        Bs=strfind(NamFileInputTemp,'/');
                        if isempty(Bs),
                            break;
                        else
                            NamFileInputTemp(1,Bs)='\';
                        end
                    end
                    
                    [NonServeS NonServeD]=system(['copy ',NamFileInputTemp, NamFileInput{jfil,2},' \\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInputTemp]);
                    
                else
                    copyfile([NamFileInput{jfil,1},NamFileInput{jfil,2}],['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInput{jfil,1}]);
                end
            end
        end
    end
end
