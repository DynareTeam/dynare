function [fOutVar,nBlockPerCPU, totCPU] = masterParallel(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar,Parallel_info,initialize)
% PARALLEL CONTEXT
% This is the most important function for the management of DYNARE parallel
% computing.
% It is the top-level function called on the master computer when parallelizing a task.

% This function have two main computational startegy for manage the matlab worker (slave process).
% 0 Simple Close/Open Stategy:
% In this case the new matlab istances (slave process) are open when
% necessary and then closed. This can happen many times during the
% simulation of a model.

% 1 Always Open Strategy:
% In this case we have a more sophisticated management of slave processes,
% which are no longer closed at the end of each job. The slave processes
% waits for a new job (if exist). If a slave do not receives a new job after a
% fixed time it is destroyed. This solution removes the computational
% time necessary to Open/Close new matlab istances.

% The first (point 0) is the default Strategy
% i.e.(Parallel_info.leaveSlaveOpen=0). This value can be changed by the
% user in xxx.mod file or it is changed by the programmer if it necessary to
% reduce the overall computational time. See for example the
% prior_posterior_statistics.m.

% The number of parallelized threads will be equal to (nBlock-fBlock+1).
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel
%  o fBlock [int]               index number of the first thread
%                               (between 1 and nBlock)
%  o nBlock [int]               index number of the last thread
%  o NamFileInput [cell array]  containins the list of input files to be
%                               copied in the working directory of remote slaves
%                               2 columns, as many lines as there are files
%                               - first column contains directory paths
%                               - second column contains filenames
%  o fname [string]             name of the function to be parallelized, and
%                               which will be run on the slaves
%  o fInputVar [struct]         structure containing local variables to be used
%                               by fName on the slaves
%  o fGlobalVar [struct]        structure containing global variables to be used
%                               by fName on the slaves
%  o Parallel_info              []
%  o initialize                 []
%
% OUTPUT
%  o fOutVar [struct vector]   result of the parallel computation, one
%                              struct per thread
%  o nBlockPerCPU [int vector] for each CPU used, indicates the number of
%                              threads run on that CPU
%  o totCPU [int]              total number of CPU used (can be lower than
%                              the number of CPU declared in "Parallel", if
%                              the number of required threads is lower)

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


% If islocal==0, create a new directory for remote computation.
% This directory is named using current data and time,
% is used only one time and then deleted.

persistent PRCDir
% PRCDir = Present Remote Computational Directory!

Strategy=Parallel_info.leaveSlaveOpen;

islocal = 1;
isHybridMatlabOctave = Parallel_info.isHybridMatlabOctave;
for j=1:length(Parallel),
    islocal=islocal*Parallel(j).Local;
end
if nargin>8 && initialize==1
    if islocal == 0,
        PRCDir=CreateTimeString();
        assignin('base','PRCDirTmp',PRCDir),
        evalin('base','options_.parallel_info.RemoteTmpFolder=PRCDirTmp;')
        evalin('base','clear PRCDirTmp,')
    else
        % Delete the traces (if existing) of last local session of computations.
        mydelete(['slaveParallel_input*.mat']);
    end
    return
end

if isfield(Parallel_info,'local_files')
    if isempty(NamFileInput),
        NamFileInput=Parallel_info.local_files;
    else
        NamFileInput=[NamFileInput;Parallel_info.local_files];
    end
end

% Deactivate some 'Parallel/Warning' message in Octave!
% Comment the line 'warning('off');' in order to view the warning message
% in Octave!

if exist('OCTAVE_VERSION'),
    warning('off');
end

% check if there are function_handles in the input or global vars when
% octave is used
if isHybridMatlabOctave || exist('OCTAVE_VERSION'),
    fInputNames = fieldnames(fInputVar);
    for j=1:length(fInputNames),
        TargetVar = fInputVar.(fInputNames{j});
        if isa(TargetVar,'function_handle'),
            TargetVar=func2str(TargetVar);
            fInputVar.(fInputNames{j})=TargetVar;
        end
    end
    
    if exist('fGlobalVar') && ~isempty(fGlobalVar),
    fInputNames = fieldnames(fGlobalVar);
    for j=1:length(fInputNames),
        TargetVar = fGlobalVar.(fInputNames{j});
        if isa(TargetVar,'function_handle'),
            TargetVar=func2str(TargetVar);
            fGlobalVar.(fInputNames{j})=TargetVar;
        end
    end
    end
end

if Strategy==1
    totCPU=0;
end


% Determine my hostname and my working directory.

DyMo=pwd;
% fInputVar.DyMo=DyMo;
if ispc, % ~(isunix || (~matlab_ver_less_than('7.4') && ismac)) ,
    [tempo, MasterName]=system('hostname');
    MasterName=deblank(MasterName);
end
% fInputVar.MasterName = MasterName;


% Save input data for use by the slaves.
switch Strategy
    case 0
        storeGlobalVars([fname,'_input.mat']);
        save([fname,'_input.mat'],'fInputVar','Parallel','-append')
        
    case 1
        if exist('fGlobalVar'),
            save(['temp_input.mat'],'fInputVar','fGlobalVar')
        else
            save(['temp_input.mat'],'fInputVar')
        end
        save(['temp_input.mat'],'Parallel','-append')
        closeSlave(Parallel,PRCDir,-1);
end


% Determine the total number of available CPUs, and the number of threads
% to run on each CPU.

[nCPU, totCPU, nBlockPerCPU, totSlaves] = distributeJobs(Parallel, fBlock, nBlock);
for j=1:totSlaves,
    PRCDirSnapshot{j}={};
end
offset0 = fBlock-1;

% Clean up remnants of previous runs.
mydelete(['comp_status_',fname,'*.mat']);
mydelete(['P_',fname,'*End.txt']);
mydelete([fname,'_output_*.mat']);
mydelete('slaveParallel_break.mat');

dynareParallelDelete([fname,'_output_*.mat'],PRCDir,Parallel);
dynareParallelDelete(['comp_status_',fname,'*.mat'],PRCDir,Parallel);
dynareParallelDelete('slaveParallel_break.mat',PRCDir,Parallel);


% Create a shell script containing the commands to launch the required
% tasks on the slaves.
fid = fopen('ConcurrentCommand1.bat','w+');


% Create the directory devoted to remote computation.
if isempty(PRCDir) && ~islocal,
    error('PRCDir not initialized!')
else
    dynareParallelMkDir(PRCDir,Parallel(1:totSlaves));
end

% Testing Zone

% 1. Display the User Strategy:

% if Strategy==0
%     disp('User Strategy Now Is Open/Close (0)');
% else
%     disp('User Strategy Now Is Always Open (1)');
% end


% 2. Display the output of 'NEW' distributeJobs.m:
%
% fBlock
% nBlock
%
%
% nCPU
% totCPU
% nBlockPerCPU
% totSlaves
%
% keyboard

% End

for j=1:totCPU,
    
    if Strategy==1
        command1 = ' ';
    end
    
    indPC=min(find(nCPU>=j));
    
    % According to the information contained in configuration file, compThread can limit MATLAB
    % to a single computational thread. By default, MATLAB makes use of the multithreading
    % capabilities of the computer on which it is running. Nevertheless
    % exsperimental results show as matlab native
    % multithreading limit the performaces when the parallel computing is active.
    
    
    if strcmp('true',Parallel(indPC).SingleCompThread),
        compThread = '-singleCompThread';
    else
        compThread = '';
    end
    
    if indPC>1
        nCPU0 = nCPU(indPC-1);
    else
        nCPU0=0;
    end
    offset = sum(nBlockPerCPU(1:j-1))+offset0;
    
    % Create a file used to monitoring if a parallel block (core)
    % computation is finished or not.
    
    fid1=fopen(['P_',fname,'_',int2str(j),'End.txt'],'w+');
    fclose(fid1);
    
    if Strategy==1,
        
        fblck = offset+1;
        nblck = sum(nBlockPerCPU(1:j));
        save temp_input.mat fblck nblck fname -append;
        copyfile('temp_input.mat',['slaveJob',int2str(j),'.mat']);
        if Parallel(indPC).Local ==0,
            fid1=fopen(['stayalive',int2str(j),'.txt'],'w+');
            fclose(fid1);
            dynareParallelSendFiles(['stayalive',int2str(j),'.txt'],PRCDir,Parallel(indPC));
            mydelete(['stayalive',int2str(j),'.txt']);
        end
        % Wait for possibly local alive CPU to start the new job or close by
        % internal criteria.
        pause(1);
        newInstance = 0;
        
        % Check if j CPU is already alive.
        if isempty(dynareParallelDir(['P_slave_',int2str(j),'End.txt'],PRCDir,Parallel(indPC)));
            fid1=fopen(['P_slave_',int2str(j),'End.txt'],'w+');
            fclose(fid1);
            if Parallel(indPC).Local==0,
                dynareParallelSendFiles(['P_slave_',int2str(j),'End.txt'],PRCDir,Parallel(indPC));
                delete(['P_slave_',int2str(j),'End.txt']);
            end
            
            newInstance = 1;
            storeGlobalVars( ['slaveParallel_input',int2str(j),'.mat']);
            save( ['slaveParallel_input',int2str(j),'.mat'],'Parallel','-append');
            % Prepare global vars for Slave.
        end
    else
        
        % If the computation is executed remotely all the necessary files
        % are created localy, then copied in remote directory and then
        % deleted (loacal)!
        
        save( ['slaveParallel_input',int2str(j),'.mat'],'j');
        
        if Parallel(indPC).Local==0,
            dynareParallelSendFiles(['P_',fname,'_',int2str(j),'End.txt'],PRCDir,Parallel(indPC));
            delete(['P_',fname,'_',int2str(j),'End.txt']);
            
            dynareParallelSendFiles(['slaveParallel_input',int2str(j),'.mat'],PRCDir,Parallel(indPC));
            delete(['slaveParallel_input',int2str(j),'.mat']);
            
        end
        
    end
    
    % DA SINTETIZZARE:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following 'switch - case' code is the core of this function!
    switch Strategy
        case 0
            
            if Parallel(indPC).Local == 1,                                  % 0.1 Run on the local machine (localhost).
                
                if ~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem), % Hybrid computing Windows <-> Unix!
                    if regexpi([Parallel(indPC).MatlabOctavePath], 'octave') % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                        command1=[Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')" &'];
                    else
                        command1=[Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')" &'];
                    end
                else    % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                    if  regexpi([Parallel(indPC).MatlabOctavePath], 'octave')
                        command1=['psexec -d -W ',DyMo, ' -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)),' -low  ',Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                    else
                        command1=['psexec -d -W ',DyMo, ' -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)),' -low  ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                    end
                end
            else                                                            % 0.2 Parallel(indPC).Local==0: Run using network on remote machine or also on local machine.
                if j==nCPU0+1,
                    dynareParallelSendFiles([fname,'_input.mat'],PRCDir,Parallel(indPC));
                    dynareParallelSendFiles(NamFileInput,PRCDir,Parallel(indPC));
                end
                
                if (~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem)), % Hybrid computing Windows <-> Unix!
                    if ispc, token='start /B ';
                    else token = '';
                    end
                    if ~isempty(Parallel(indPC).Port),
                        ssh_token = ['-p ',Parallel(indPC).Port];
                    else
                        ssh_token = '';
                    end
                    % To manage the diferences in Unix/Windows OS syntax.
                    remoteFile=['remoteDynare',int2str(j)];
                    fidRemote=fopen([remoteFile,'.m'],'w+');
                    if regexpi([Parallel(indPC).MatlabOctavePath], 'octave'),% Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                        remoteString=['default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
                        command1=[token, 'ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' "cd ',Parallel(indPC).RemoteDirectory,'/',PRCDir, '; ',Parallel(indPC).MatlabOctavePath,' -f --eval ',remoteFile,' " &'];
                    else
                        remoteString=['addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
                        command1=[token, 'ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' "cd ',Parallel(indPC).RemoteDirectory,'/',PRCDir, '; ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r ',remoteFile,';" &'];
                    end
                    fprintf(fidRemote,'%s\n',remoteString);
                    fclose(fidRemote);
                    dynareParallelSendFiles([remoteFile,'.m'],PRCDir,Parallel(indPC));
                    delete([remoteFile,'.m']);
                else
                    if ~strcmpi(Parallel(indPC).ComputerName,MasterName),  % 0.3 Run on a remote machine!
                        % Hybrid computing Matlab(Master)-> Octave(Slaves) and Vice Versa!
                        if  regexpi([Parallel(indPC).MatlabOctavePath], 'octave')
                            command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -u ',Parallel(indPC).UserName,' -p ',Parallel(indPC).Password,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                ' -low  ',Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                        else
             
                            command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -u ',Parallel(indPC).UserName,' -p ',Parallel(indPC).Password,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                ' -low  ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                        end
                    else                                                  % 0.4 Run on the local machine via the network
                        % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                        if  regexpi([Parallel(indPC).MatlabOctavePath], 'octave')
                            command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                ' -low  ',Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                        else
                            command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                ' -low  ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                        end
                    end
                end
            end
            
            
        case 1
            if Parallel(indPC).Local == 1 && newInstance,                       % 1.1 Run on the local machine.
                if (~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem)),  % Hybrid computing Windows <-> Unix!                   
                    if regexpi([Parallel(indPC).MatlabOctavePath], 'octave')    % Hybrid computing Matlab(Master)-> Octave(Slaves) and Vice Versa!
                        command1=[Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')" &'];
                    else
                        command1=[Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')" &'];
                    end
                else    % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                    if  regexpi([Parallel(indPC).MatlabOctavePath], 'octave')
                        command1=['psexec -d -W ',DyMo, ' -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)),' -low  ',Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7'');addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                    else
                        command1=['psexec -d -W ',DyMo, ' -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)),' -low  ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                    end
                end
            elseif Parallel(indPC).Local==0,                                % 1.2 Run using network on remote machine or also on local machine.
                if j==nCPU0+1,
                    dynareParallelSendFiles(NamFileInput,PRCDir,Parallel(indPC));
                end
                dynareParallelSendFiles(['P_',fname,'_',int2str(j),'End.txt'],PRCDir,Parallel(indPC));
                delete(['P_',fname,'_',int2str(j),'End.txt']);
                if newInstance,
                    dynareParallelSendFiles(['slaveJob',int2str(j),'.mat'],PRCDir,Parallel(indPC));
                    delete(['slaveJob',int2str(j),'.mat']);
                    dynareParallelSendFiles(['slaveParallel_input',int2str(j),'.mat'],PRCDir,Parallel(indPC))
                    if (~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem)), % Hybrid computing Windows <-> Unix!
                        if ispc, token='start /B ';
                        else token = '';
                        end
                        if ~isempty(Parallel(indPC).Port),
                            ssh_token = ['-p ',Parallel(indPC).Port];
                        else
                            ssh_token = '';
                        end
                        % To manage the diferences in Unix/Windows OS syntax.
                        remoteFile=['remoteDynare',int2str(j)];
                        fidRemote=fopen([remoteFile,'.m'],'w+');
                        if regexpi([Parallel(indPC).MatlabOctavePath], 'octave') % Hybrid computing Matlab(Master)-> Octave(Slaves) and Vice Versa!
                            remoteString=['default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),');'];
                            command1=[token, 'ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' "cd ',Parallel(indPC).RemoteDirectory,'/',PRCDir '; ',Parallel(indPC).MatlabOctavePath,' -f --eval ',remoteFile,' " &'];
                        else
                            remoteString=['addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),');'];
                            command1=[token, 'ssh ',ssh_token,' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' "cd ',Parallel(indPC).RemoteDirectory,'/',PRCDir '; ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r ',remoteFile,';" &'];
                        end
                        fprintf(fidRemote,'%s\n',remoteString);
                        fclose(fidRemote);
                        dynareParallelSendFiles([remoteFile,'.m'],PRCDir,Parallel(indPC));
                        delete([remoteFile,'.m']);
                    else
                        if ~strcmpi(Parallel(indPC).ComputerName,MasterName), % 1.3 Run on a remote machine.
                            % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                            if  regexpi([Parallel(indPC).MatlabOctavePath], 'octave')
                                command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -u ',Parallel(indPC).UserName,' -p ',Parallel(indPC).Password,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                    ' -low  ',Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7'');addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                            else
                                command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -u ',Parallel(indPC).UserName,' -p ',Parallel(indPC).Password,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                    ' -low  ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                            end
                        else                                                % 1.4 Run on the local machine via the network.
                            % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                            if  regexpi([Parallel(indPC).MatlabOctavePath], 'octave')
                                command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                    ' -low  ',Parallel(indPC).MatlabOctavePath,' -f --eval "default_save_options(''-v7''); addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                            else
                                command1=['psexec \\',Parallel(indPC).ComputerName,' -d -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\ -a ',int2str(Parallel(indPC).CPUnbr(j-nCPU0)), ...
                                    ' -low  ',Parallel(indPC).MatlabOctavePath,' -nosplash -nodesktop -minimize ',compThread,' -r "addpath(''',Parallel(indPC).DynarePath,'''), dynareroot = dynare_config(); slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                            end
                        end
                    end
                else
                    % When the user user strategy is equal to 1, you must
                    % do PRCDirSnapshot here to to avoid problems of
                    % synchronization.
                    
                    if isempty(PRCDirSnapshot{indPC}),
                        PRCDirSnapshot(indPC)=dynareParallelSnapshot(PRCDir,Parallel(indPC));
                        PRCDirSnapshotInit(indPC) = PRCDirSnapshot(indPC);
                    else
                        PRCDirSnapshot(indPC)=dynareParallelGetNewFiles(PRCDir,Parallel(indPC),PRCDirSnapshot(indPC));
                    end
                    dynareParallelSendFiles(['slaveJob',int2str(j),'.mat'],PRCDir,Parallel(indPC));
                    delete(['slaveJob',int2str(j),'.mat']);
                    
                end
            end
            
    end
    
    fprintf(fid,'%s\n',command1);
    
end

% In This way we are sure that the file 'ConcurrentCommand1.bat' is
% closed and then it can be deleted!
while (1)
    StatusOfCC1_bat = fclose(fid);
    if StatusOfCC1_bat==0
        break
    end
end
% Snapshot  of the contents of all the directories involved in parallel
% computing. This is necessary when I want to copy continuously the files produced by
% the slaves ...
% If the compuation is 'Local' it is not necessary to do it ...

if Strategy==0 || newInstance, % See above.
    PRCDirSnapshot=dynareParallelSnapshot(PRCDir,Parallel(1:totSlaves));
    PRCDirSnapshotInit = PRCDirSnapshot;
    
    % Run the slaves.
    if  ~ispc, %isunix || (~matlab_ver_less_than('7.4') && ismac),
        system('sh ConcurrentCommand1.bat &');
        pause(1)
    else
        
        if exist('OCTAVE_VERSION')
            % Redirect the standard output to the file 'OctaveStandardOutputMessage.txt'!
            % This file is saved in the Model directory.
            system('ConcurrentCommand1.bat > OctaveStandardOutputMessage.txt');
        else
            system('ConcurrentCommand1.bat');
        end
    end
end


% For matlab enviroment with options_.console_mode = 0:
% create a parallel (local/remote) specialized computational status bars!

global options_


% Create a parallel (local/remote) specialized computational status bars!

if exist('OCTAVE_VERSION') || (options_.console_mode == 1),
    diary off;
    if exist('OCTAVE_VERSION')
        printf('\n');
    else
        fprintf('\n');
    end
else
    hfigstatus = figure('name',['Parallel ',fname],...
        'DockControls','off', ...
        'IntegerHandle','off', ...
        'Interruptible','off', ...
        'MenuBar', 'none', ...
        'NumberTitle','off', ...
        'Renderer','Painters', ...
        'Resize','off');
    
    vspace = 0.1;
    ncol = ceil(totCPU/10);
    hspace = 0.9/ncol;
    hstatus(1) = axes('position',[0.05/ncol 0.92 0.9/ncol 0.03], ...
        'box','on','xtick',[],'ytick',[],'xlim',[0 1],'ylim',[0 1]);
    set(hstatus(1),'Units','pixels')
    hpixel = get(hstatus(1),'Position');
    hfigure = get(hfigstatus,'Position');
    hfigure(4)=hpixel(4)*10/3*min(10,totCPU);
    set(hfigstatus,'Position',hfigure)
    set(hstatus(1),'Units','normalized'),
    vspace = max(0.1,1/totCPU);
    vstart = 1-vspace+0.2*vspace;
    for j=1:totCPU,
        jrow = mod(j-1,10)+1;
        jcol = ceil(j/10);
        hstatus(j) = axes('position',[0.05/ncol+(jcol-1)/ncol vstart-vspace*(jrow-1) 0.9/ncol 0.3*vspace], ...
            'box','on','xtick',[],'ytick',[],'xlim',[0 1],'ylim',[0 1]);
        hpat(j) = patch([0 0 0 0],[0 1 1 0],'r','EdgeColor','r');
        htit(j) = title(['Initialize ...']);
        
    end
    
    cumBlockPerCPU = cumsum(nBlockPerCPU);
end
pcerdone = NaN(1,totCPU);
idCPU = NaN(1,totCPU);



% Wait for the slaves to finish their job, and display some progress
% information meanwhile.

% Caption for console mode computing ...

if (options_.console_mode == 1) ||  exist('OCTAVE_VERSION')
    
    if ~exist('OCTAVE_VERSION')
        if strcmpi([Parallel(indPC).MatlabOctavePath], 'octave')
            RjInformation='Hybrid Computing Is Active: Remote jobs are computed by Octave!';
            fprintf([RjInformation,'\n\n']);
        end
    end
    
    fnameTemp=fname;
    
    L=length(fnameTemp);
    
    PoCo=strfind(fnameTemp,'_core');
    
    for i=PoCo:L
        if i==PoCo
            fnameTemp(i)=' ';
        else
            fnameTemp(i)='.';
        end
    end
    
    for i=1:L
        if  fnameTemp(i)=='_';
            fnameTemp(i)=' ';
        end
    end
    
    fnameTemp(L)='';
    
    Information=['Parallel ' fnameTemp ' Computing ...'];
    if exist('OCTAVE_VERSION')
        if (~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem)) && (Strategy==0)
            printf('\n');
            pause(2);
        end
        
        printf([Information,'\n\n']);
    else
        fprintf([Information,'\n\n']);
    end
    
end


% Testing Zone

% Check the new copy file strategy ...
global NuoviFilecopiati
NuoviFilecopiati=zeros(1,totSlaves);
% End

ForEver=1;
statusString = '';
flag_CloseAllSlaves=0;

while (ForEver)
    
    waitbarString = '';
    statusString0 = repmat('\b',1,length(sprintf(statusString, 100 .* pcerdone)));
    statusString = '';
    
    pause(1)
    
    try
        if islocal ==0,
            dynareParallelGetFiles(['comp_status_',fname,'*.mat'],PRCDir,Parallel(1:totSlaves));
        end
    catch
    end
    
    for j=1:totCPU,
        try
            if ~isempty(['comp_status_',fname,int2str(j),'.mat'])
                load(['comp_status_',fname,int2str(j),'.mat']);
%                 whoCloseAllSlaves = who(['comp_status_',fname,int2str(j),'.mat','CloseAllSlaves']);
                if exist('CloseAllSlaves') && flag_CloseAllSlaves==0,
                    flag_CloseAllSlaves=1;
                    whoiamCloseAllSlaves=j;
                    closeSlave(Parallel(1:totSlaves),PRCDir,1);
                end
            end
            pcerdone(j) = prtfrc;
            idCPU(j) = njob;
            if exist('OCTAVE_VERSION') || (options_.console_mode == 1),
                if (~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem))
                    statusString = [statusString, int2str(j), ' %3.f%% done! '];
                else
                    statusString = [statusString, int2str(j), ' %3.f%% done! '];
                end
            else
                status_String{j} = waitbarString;
                status_Title{j} = waitbarTitle;
            end
        catch % ME
            % To define!
            if exist('OCTAVE_VERSION') || (options_.console_mode == 1),
                if (~ispc || strcmpi('unix',Parallel(indPC).OperatingSystem))
                    statusString = [statusString, int2str(j), ' %3.f%% done! '];
                else
                    statusString = [statusString, int2str(j), ' %3.f%% done! '];
                end
            end
        end
    end
    if exist('OCTAVE_VERSION') || (options_.console_mode == 1),
        if exist('OCTAVE_VERSION')
            printf([statusString,'\r'], 100 .* pcerdone);
        else
            if ~isempty(statusString)
                fprintf([statusString0,statusString], 100 .* pcerdone);
            end
        end
        
    else
        for j=1:totCPU,
            try
                set(hpat(j),'XData',[0 0 pcerdone(j) pcerdone(j)]);
                set(htit(j),'String',[status_Title{j},' - ',status_String{j}]);
            catch
                
            end
        end
    end
    
    % Check if the slave(s) has generated some new files remotely.
    % 1. The files .log and .txt are not copied.
    % 2. The comp_status_*.mat files are managed separately.
    
    if exist('OCTAVE_VERSION'), % to avoid synchronism problems
        try
            PRCDirSnapshot=dynareParallelGetNewFiles(PRCDir,Parallel(1:totSlaves),PRCDirSnapshot);
        catch
        end
    else
        PRCDirSnapshot=dynareParallelGetNewFiles(PRCDir,Parallel(1:totSlaves),PRCDirSnapshot);
    end
    
    if isempty(dynareParallelDir(['P_',fname,'_*End.txt'],PRCDir,Parallel(1:totSlaves)));
        HoTuttiGliOutput=0;
        for j=1:totCPU,
            
            % Checking if the remote computation is finished and if we copied all the output here.
            if ~isempty(dir([fname,'_output_',int2str(j),'.mat']))
                HoTuttiGliOutput=HoTuttiGliOutput+1;
            else
                indPC=min(find(nCPU>=j));
                dynareParallelGetFiles([fname,'_output_',int2str(j),'.mat'],PRCDir,Parallel(indPC));
            end
        end
        
        if HoTuttiGliOutput==totCPU,
            mydelete(['comp_status_',fname,'*.mat']);
            if exist('OCTAVE_VERSION')|| (options_.console_mode == 1),
                if exist('OCTAVE_VERSION')
                    printf('\n');
                    printf(['End Parallel Session ....','\n\n']);
                else
                    fprintf('\n');
                    fprintf(['End Parallel Session ....','\n\n']);
                end
                diary on;
            else
                close(hfigstatus),
            end
            
            break
        else
            disp('Waiting for output files from slaves ...')
        end
    end
    
end


% Load and format remote output.
iscrash = 0;
PRCDirSnapshot=dynareParallelGetNewFiles(PRCDir,Parallel(1:totSlaves),PRCDirSnapshot);

for j=1:totCPU,
    indPC=min(find(nCPU>=j));
    load([fname,'_output_',int2str(j),'.mat'],'fOutputVar');
    delete([fname,'_output_',int2str(j),'.mat']);
    if isfield(fOutputVar,'OutputFileName') && Parallel(indPC).Local==0,
        %   Check if input files have been updated!
        OutputFileName=fOutputVar.OutputFileName;        
        tmp0='';
        for i=1:size(NamFileInput,1),
            FileList = regexp(strrep(PRCDirSnapshot{indPC},'\','/'),strrep(strrep([NamFileInput{i,:}],'\','/'),'*','(\w*)'),'match');
            for k=1:length(FileList),
                if ~isempty(FileList{k}),
                    if isempty(tmp0),
                        tmp0=FileList{k}{1};
                    else
                        tmp0=char(tmp0,FileList{k}{1});
                    end
                end
            end
        end
        for i=1:size(OutputFileName,1),
            tmp1='';
            FileList = regexp(cellstr(tmp0),strrep(strrep([OutputFileName{i,:}],'\','/'),'*','(\w*)'),'match');
            FileList0 = regexp(cellstr(tmp0),strrep([OutputFileName{i,2}],'*','(\w*)'),'match');
            for k=1:length(FileList),
                if ~isempty(FileList{k}),
                    if isempty(tmp1),
                        tmp1=FileList0{k}{1};
                    else
                        tmp1=char(tmp1,FileList0{k}{1});
                    end
                end
            end
            for k=1:size(tmp1,1),
                    dynareParallelGetFiles([OutputFileName(i,1) {tmp1(k,:)}],PRCDir,Parallel(indPC));
            end
        end
        % check if some output file is missing!
        for i=1:size(OutputFileName,1),
            tmp1=dynareParallelDir([OutputFileName{i,:}],PRCDir,Parallel(indPC));
            tmp1 = regexp(cellstr(tmp1),strrep([OutputFileName{i,2}],'*','(\w*)'),'match');
            tmp1 = char(tmp1{:});
            tmp2=ls([OutputFileName{i,:}]);
            for ij=1:size(tmp1,1),
                icheck = regexp(cellstr(tmp2),tmp1(ij,:),'once');
                isOutputFileMissing=1;
                for ik=1:size(tmp2,1),
                    if ~isempty(icheck{ik}),
                        isOutputFileMissing=0;
                    end
                end
                if isOutputFileMissing,
                    dynareParallelGetFiles([OutputFileName(i,1) {tmp1(ij,:)}],PRCDir,Parallel(indPC));
                end
            end
        end

    end
    if isfield(fOutputVar,'error'),
        disp(['Job number ',int2str(j),' crashed with error:']);
        iscrash=1;
        disp([fOutputVar.error.message]);
        for jstack=1:length(fOutputVar.error.stack)
            fOutputVar.error.stack(jstack),
        end
    elseif flag_CloseAllSlaves==0,
        fOutVar(j)=fOutputVar;
    elseif j==whoiamCloseAllSlaves,
        fOutVar=fOutputVar;        
    end
end

if flag_CloseAllSlaves==1,
    closeSlave(Parallel(1:totSlaves),PRCDir,-1);
end

if iscrash,
    error('Remote jobs crashed');
end

pause(1), % Wait for all remote diary off completed

% Cleanup.
dynareParallelGetFiles('*.log',PRCDir,Parallel(1:totSlaves));

switch Strategy
    case 0
        for indPC=1:min(find(nCPU>=totCPU)),
            if Parallel(indPC).Local == 0
                dynareParallelRmDir(PRCDir,Parallel(indPC));
            end
            
            if isempty(dir('dynareParallelLogFiles'))
                [A B C]=rmdir('dynareParallelLogFiles');
                mkdir('dynareParallelLogFiles');
            end
            try
                copyfile('*.log','dynareParallelLogFiles');
                mydelete([fname,'*.log']);
            catch
            end
            
            mydelete(['*_core*_input*.mat']);
            
        end
        
        delete ConcurrentCommand1.bat
    case 1
        delete(['temp_input.mat'])
        if newInstance,
            if isempty(dir('dynareParallelLogFiles'))
                [A B C]=rmdir('dynareParallelLogFiles');
                mkdir('dynareParallelLogFiles');
            end
        end
        copyfile('*.log','dynareParallelLogFiles');
        if newInstance,
            delete ConcurrentCommand1.bat
        end
        dynareParallelDelete(['comp_status_',fname,'*.mat'],PRCDir,Parallel);
        for indPC=1:min(find(nCPU>=totCPU)),
            if Parallel(indPC).Local == 0,
                dynareParallelDeleteNewFiles(PRCDir,Parallel(indPC),PRCDirSnapshotInit(indPC),'.log');
                for ifil=1:size(NamFileInput,1),
                    dynareParallelDelete([NamFileInput{ifil,:}],PRCDir,Parallel(indPC));
                end
            end
        end
end




