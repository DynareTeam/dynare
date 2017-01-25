function [nCPU, totCPU, nBlockPerCPU, totSLAVES] = distributeJobs(Parallel, fBlock, nBlock)
% PARALLEL CONTEXT
% In parallel context this function is used to determine the total number of available CPUs,
% and the number of threads to run on each CPU.
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel
%  o fBlock [int]               index number of the first job (e.g. MC iteration or MH block)
%                               (between 1 and nBlock)
%  o nBlock [int]               index number of the last job.
%
% OUTPUT
%  o nBlockPerCPU [int vector]  for each CPU used, indicates the number of
%                               threads run on that CPU
%  o totCPU [int]               total number of CPU used (can be lower than
%                               the number of CPU declared in "Parallel", if
%                               the number of required threads is lower!)
%  o nCPU                       the number of CPU in user format.
%  o totSLAVES                  the number of cluster's node currently
%                               involved in parallel computing step.
%                               It is a number between 1 and length(Parallel).


% Copyright (C) 2010-2011 Dynare Team
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


% The Parallel vector has already been sorted
% (in accord with the CPUWeight values) in DESCENDING order in 
% InitializeComputationalEnvironment!

totCPU=0;

lP=length(Parallel);
CPUWeight=ones(1,length(Parallel))*(-1);

for j=1:lP,    
    if mod(length(Parallel(j).CPUnbr),Parallel(j).NumberOfThreadsPerJob)
        skipline()
        disp(['PARALLEL_ERROR:: NumberOfThreadsPerJob = ',int2str(Parallel(j).NumberOfThreadsPerJob),' is not an exact divisor of number of CPUs = ',int2str(length(Parallel(j).CPUnbr)),'!'])
        disp(['                 You must re-set properly NumberOfThreadsPerJob of node ' int2str(j) ' ' Parallel(j).ComputerName])
        disp(['                 in your configuration file'])
        error(['PARALLEL_ERROR:: NumberOfThreadsPerJob is not an exact divisor of CPUnbr'])
    end
    nCPU(j)=length(Parallel(j).CPUnbr)/Parallel(j).NumberOfThreadsPerJob;
    totCPU=totCPU+nCPU(j);    
    CPUWeight(j)=str2num(Parallel(j).NodeWeight);
end


% Copy of original nCPU.
nCPUoriginal=nCPU;

nCPU=cumsum(nCPU);


% Number of Nodes in Cluster.
nC=lP;

% Numbers of Jobs.
NumbersOfJobs=nBlock-fBlock+1;

SumOfJobs=0;
JobsForNode=zeros(1,nC);

for j=1:lP,
    CPUWeight(j)=str2num(Parallel(j).NodeWeight)*nCPUoriginal(j);
end
CPUWeight=CPUWeight./sum(CPUWeight);

% Redistributing the jobs among the cluster nodes according to the
% CPUWeight.
for i=1:nC
    
    JobsForNode(i)=CPUWeight(i)*NumbersOfJobs;
    
    % Many choices are possible:
    
    % JobsForNode(i)=round(JobsForNode(i));
    % JobsForNode(i)=floor(JobsForNode(i));
      JobsForNode(i)=ceil(JobsForNode(i));
    
end

% Check if there are more (fewer) jobs.
% This can happen when we use ceil, round, ... functions.
SumOfJobs=sum(JobsForNode);

if SumOfJobs~=NumbersOfJobs
    
    if SumOfJobs>NumbersOfJobs
        
        % Many choices are possible:
        
        % - Remove the excess works at the node that has the greatest
        %   number of jobs.
        % - Remove the excess works at the node slower.
        
        VerySlow=nC;
        
        while SumOfJobs>NumbersOfJobs
            
            if JobsForNode(VerySlow)==0
                VerySlow=VerySlow-1;
                continue
            end
            JobsForNode(VerySlow)=JobsForNode(VerySlow)-1;
            SumOfJobs=SumOfJobs-1;
        end
        
    end
    
    if SumOfJobs<NumbersOfJobs
        
        % Many choices are possible:
        % - ... (see above).
        
        [NonServe VeryFast]= min(CPUWeight);
        
        while SumOfJobs<NumbersOfJobs
            JobsForNode(VeryFast)=JobsForNode(VeryFast)+1;
            SumOfJobs=SumOfJobs+1;
        end
        
    end
end


% Redistributing the jobs among the cpu/core nodes.

JobsForCpu=zeros(1,nCPU(nC));
JobAssignedCpu=0;

RelativePosition=1;

for i=1:nC
    
    % Many choices are possible:
    % - ... (see above).
     
    JobAssignedCpu=max(1,floor(JobsForNode(i)/nCPUoriginal(i)));
    
    ChekOverFlow=0;
    
    for j=RelativePosition:nCPU(i)
        JobsForCpu(j)=JobAssignedCpu;
        ChekOverFlow=ChekOverFlow+JobAssignedCpu;
        
        if ChekOverFlow>=JobsForNode(i)
            break;
        end
        
    end
    
    % Check if there are more (fewer) jobs.
    % This can happen when we use ceil, round, ... functions.
    
    if ChekOverFlow ~=(JobsForNode(i))
        
        if ChekOverFlow >(JobsForNode(i))
            while ChekOverFlow>JobsForNode(i)
                JobsForCpu(nCPU(i))=JobsForCpu(nCPU(i))-1;
                ChekOverFlow=ChekOverFlow-1;
            end
        end
        
        if ChekOverFlow <(JobsForNode(i))
            while ChekOverFlow<JobsForNode(i)
                JobsForCpu(nCPU(i))=JobsForCpu(nCPU(i))+1;
                ChekOverFlow=ChekOverFlow+1;
            end
        end
    end
    
    RelativePosition=nCPU(i)+1;
    
end

% Reshape the variables totCPU,totSLAVES and nBlockPerCPU in accord with
% the syntax rquired by a previous version of parallel package ...

totCPU=0;
totSLAVES=0;
nBlockPerCPU=[];

for i=1:nCPU(nC)
    if JobsForCpu(i)~=0
        totCPU=totCPU+1;
    end
end

for i=1:nC
    if JobsForNode(i)~=0;
        totSLAVES=totSLAVES+1;
    end
end

RelativeCounter=1;
for i=1:nCPU(nC)
    if JobsForCpu(i)~=0
        nBlockPerCPU(RelativeCounter)=JobsForCpu(i);
        RelativeCounter=RelativeCounter+1;
    end
end


