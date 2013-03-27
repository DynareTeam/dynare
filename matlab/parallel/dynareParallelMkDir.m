function dynareParallelMkDir(PRCDir,Parallel)
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



if nargin ==0,
    disp('dynareParallelMkDir(dirname,Parallel)')
    return
end

for indPC=1:length(Parallel)
    if Parallel(indPC).Local==0,
        if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem), %isunix || (~matlab_ver_less_than('7.4') && ismac),
            if ~isempty(Parallel(indPC).Port),
                ssh_token = ['-p ',Parallel(indPC).Port];
            else
                ssh_token = '';
            end
            [NonServeS NonServeD]=system(['ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' mkdir -p ',Parallel(indPC).RemoteDirectory,'/',PRCDir]);
        else
            [NonServeS NonServeD]=mkdir(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir]);
        end
    end
end

return