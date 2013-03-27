function dynareParallelRmDir(PRCDir,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized version of rmdir() function.
%
% INPUTS
%  o PRCDir         []   ... 
%  o Parallel       []   ...  
%
%  OUTPUTS
%  None
%
%
%
% Copyright (C) 2009-2012 Dynare Team
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



if nargin ==0,
    disp('dynareParallelRmDir(fname)')
    return
end

% security check of remote folder delete
ok(1)=isempty(strfind(PRCDir,'..'));
tmp1=strfind(PRCDir,'2');
ok(2)=tmp1(1)==1;
ok(3)=~isempty(strfind(PRCDir,'-'));
ok(4)=~isempty(strfind(PRCDir,'h'));
ok(5)=~isempty(strfind(PRCDir,'m'));
ok(6)=~isempty(strfind(PRCDir,'s'));
ok(7)=~isempty(PRCDir);

if sum(ok)<7,
    error('The name of the remote tmp folder does not comply the security standards!'),
end

for indPC=1:length(Parallel),
    ok(1)=isempty(strfind(Parallel(indPC).RemoteDirectory,'..'));
    if sum(ok)<7,
        error('The remote folder path structure does not comply the security standards!'),
    end
    while (1)
        if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem),
            if ~isempty(Parallel(indPC).Port),
                ssh_token = ['-p ',Parallel(indPC).Port];
            else
                ssh_token = '';
            end
            [stat NonServe] = system(['ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' rm -fr ',Parallel(indPC).RemoteDirectory,'/',PRCDir,]);
            break;
        else
            if exist('OCTAVE_VERSION'), % Patch for peculiar behaviour of rmdir under Windows.
                                        % It is necessary because the command rmdir always ask at the user to confirm your decision before
                                        % deleting a directory: this stops the computation! The Octave native function 'confirm_recursive_rmdir'
                                        % disable this mechanism.
                val = confirm_recursive_rmdir (false);
                [stat, mess, id] = rmdir(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir],'s');
                
            else
                [stat, mess, id] = rmdir(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir],'s');
            end
            
            if stat==1,
                break,
            else
                if isempty(dynareParallelDir(PRCDir,'',Parallel));
                    break,
                else
                    pause(1);
                end
            end
        end
    end
end

return