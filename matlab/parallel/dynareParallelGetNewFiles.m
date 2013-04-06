function [PRCDirSnapshot]=dynareParallelGetNewFiles(PRCDir,Parallel,PRCDirSnapshot)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized function able to ...
%
%
% INPUTS
%
%  o PRCDir           []   ...
%  o Parallel         []   ...
%  o PRCDirSnapshot   []   ...
%
%
%  OUTPUTS
%  o PRCDirSnapshot   []   ...
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


NewFilesFromSlaves={};

% try
for indPC=1:length(Parallel),
    
    if Parallel(indPC).Local==0;
        [NewFilesFromSlaves, PRCDirSnapshot{indPC}]=dynareParallelFindNewFiles(PRCDirSnapshot{indPC},Parallel(indPC), PRCDir);
        if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem),
            fS='/';
        else
            fS='\';
        end

        if ~isempty(NewFilesFromSlaves)
            pause(1);

            for i=1:length(NewFilesFromSlaves)
                SlashNumberAndPosition=[];
                PRCDirPosition=findstr(NewFilesFromSlaves{i}, ([PRCDir]));
                sT=NewFilesFromSlaves{i};
                sT(1:(PRCDirPosition+length([PRCDir]))-2)=[];
                sT(1)='.';
                SlashNumberAndPosition=findstr(sT,fS);
                fileaddress={sT(1:SlashNumberAndPosition(end)),sT(SlashNumberAndPosition(end)+1:end)};
                dynareParallelGetFiles(fileaddress,PRCDir,Parallel(indPC));

                disp('New file copied in local -->');
                disp(fileaddress{2});
                disp('<--');

            end
        else
            continue
        end


    end
end



