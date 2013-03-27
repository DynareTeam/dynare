function closeSlave(Parallel,TmpFolder,partial),
% PARALLEL CONTEXT
% In parallel context, this utility closes all remote matlab instances
% called by masterParallel when strategy (1) is active i.e. always open (which leaves
% open remote matlab instances).
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel.
%  o TmpFolder        string    if islocal==0, is the name of didectory devoted to remote computation.
%                               This directory is named using current date
%                               and is used only one time and then deleted.
%                               If islocal==1, TmpFolder=''.
% 
%
% OUTPUTS
%   None
%
% Copyright (C) 2010-2012 Dynare Team
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

if nargin<3,
    partial=0;
end

s=warning('off');

if partial==1
    save('slaveParallel_break.mat','partial')
    for indPC=1:length(Parallel),
        if (Parallel(indPC).Local==0),
            dynareParallelSendFiles('slaveParallel_break.mat',TmpFolder,Parallel(indPC));
        end
    end
%     delete('slaveParallel_break')
    return
end
if partial==-1
    delete('slaveParallel_break.mat')
    for indPC=1:length(Parallel),
        if (Parallel(indPC).Local==0),
            dynareParallelDelete( 'slaveParallel_break.mat',TmpFolder,Parallel(indPC));
        end
    end
%     delete('slaveParallel_break')
    return
end

for indPC=1:length(Parallel),
    if (Parallel(indPC).Local==0),
        dynareParallelDelete( 'slaveParallel_input*.mat',TmpFolder,Parallel(indPC));
    end
    
   
    delete( 'slaveParallel_input*.mat');
    delete( 'slaveJob*.mat');
    pause(1)
    delete(['slaveParallel_*.log']);
    delete ConcurrentCommand1.bat;
    
end

while(1)
    if isempty(dynareParallelDir(['P_slave_',int2str(j),'End.txt'],TmpFolder,Parallel));
        for indPC=1:length(Parallel),
            if (Parallel(indPC).Local==0),
                dynareParallelRmDir(TmpFolder,Parallel(indPC)),
            end
        end
        break
        
    end
end

s=warning('on');
 
