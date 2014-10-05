function [fOutVar,nBlockPerCPU, CPUsinUse] = masterParallel2(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar,Parallel_info,initialize)

%% pack global variables, send to core function via input arguments

 globalVars = who('global');  
 
 for j=1:length(globalVars)
        eval(['global ',globalVars{j},';'])    
 end 
       
 for i=1:length(globalVars)
     name=globalVars{i};
     globalVars2.(globalVars{i})=eval(name);
 end
  
 fInputVar.global=globalVars2;
 workdir=cd;
 fInputVar.HostDir=workdir;
 HostNr=find([(options_.parallel.Local)]);
 Parallel(HostNr).ComputerName=char(java.net.InetAddress.getLocalHost.getHostName);
 fInputVar.Parallel = Parallel;
 
%% retrieve cluster information from config file

 totCPU=0;
 lP=length(Parallel);

 for j=1:lP,
    nCPU(j)=length(Parallel(j).CPUnbr);
    totCPU=totCPU+nCPU(j);
 end

 NumberOfBlocks=nBlock-fBlock+1;
 remainder=0;
 
 switch rem(NumberOfBlocks,totCPU)
     
     case 0                                            % even distribution across CPUs
          BlocksPerCPU=NumberOfBlocks/totCPU;
          CPUsinUse=totCPU;
          
     case NumberOfBlocks                               % more CPUs than blocks, each block on a different CPU
          BlocksPerCPU=1;
          CPUsinUse=NumberOfBlocks;
     
     otherwise
          BlocksPerCPU=floor(NumberOfBlocks/totCPU);
          CPUsinUse=totCPU;
          remainder=rem(NumberOfBlocks,totCPU);        % the first |rem| workers get an extra chain
 end
 
 %% slice the task for workers
 slices=zeros(CPUsinUse,2);
 nBlockPerCPU=zeros(1,CPUsinUse);
 fblck=fBlock;

 for i=1:CPUsinUse
    nblck=(fblck+BlocksPerCPU)-1;
    if remainder && ( i <= remainder)                      % the first |rem| workers get an extra chain
       nblck=nblck+1;                 
    end
    slices(i,1)=fblck;
    slices(i,2)=nblck;
    nBlockPerCPU(1,i)=(nblck-fblck)+1;
    fblck=nblck+1;
 end
 
%% 
if ~isdir([ workdir '/diary/'])
    mkdir([ workdir '/diary/'])
end

%%
if ( (options_.Cluster_settings == 1) || (options_.Cluster_settings == 3) ) % set up cluster environment, maximum nuber of CPUs as defined in config-file
    
    if verLessThan('matlab', '8')
        cluster = findResource('scheduler', 'type', 'Torque');
        cluster.ClusterMatlabRoot = Parallel(1).MatlabOctavePath;
        cluster.ClusterSize = totCPU;
        cluster.ResourceTemplate = '-l nodes=^N^';
        cluster.DataLocation = [ workdir '/diary/'];
    else
        cluster = parallel.cluster.Torque();
        cluster.ClusterMatlabRoot = '/opt/matlab2012b';                    
        cluster.JobStorageLocation =[ workdir '/diary/'];
        cluster.NumWorkers = totCPU;
    end
    
    if options_.Cluster_settings == 1
        cluster.HasSharedFilesystem = true;
    elseif options_.Cluster_settings == 3% store temporarily on node
        cluster.HasSharedFilesystem = false;
    end
    
    cluster.RcpCommand = 'scp';
    cluster.RshCommand = 'ssh';
    cluster.SubmitArguments = '-q verylong -l walltime=12:00:00 -m bea -M youremail'; % enter email
    
elseif options_.Cluster_settings == 2 % set up worker environment -- local
    
    parallel.defaultClusterProfile('local');
    cluster = parcluster();
    cluster.JobStorageLocation = [ workdir '/diary/'];
    
else
    error('Unknown parallel option')
end

 %% attach files to input-arguments, if no shared filesystem available
 if (options_.Cluster_settings == 3) 
    SendFiles={};
    counter=1;
    for i=1:size(NamFileInput,1)
        files=dir([NamFileInput{i,1} NamFileInput{i,2}]);
        for j=1:size(files,1)
            fid=fopen([NamFileInput{i,1} files(j,1).name]);
            SendFiles{counter}.data=fread(fid);
            fclose(fid);
            SendFiles{counter}.name=[files(j,1).name];
            counter=counter+1;
        end
    end
    fInputVar.SendFiles=SendFiles;
 end
 
 %% job-dispatch
 
 job = createJob(cluster);
 
  pathfolders=regexp(path,':','split'); % current path settings of localhost
  workfolders=regexp(genpath(cd),':','split');
  combined={pathfolders{:}, workfolders{:}};
  set(job,'AdditionalPaths',combined)
  
 for i=1:CPUsinUse
    
     a=createTask(job, fname, 1, {fInputVar,slices(i,1),slices(i,2),i});  
     a.CaptureDiary = true;
     
 end
 
 submit(job);
 display(CPUsinUse)
 display('(Number of submitted jobs)')
 display(NumberOfBlocks)
 display(slices)
 get(job,'Tasks')
 state=get(job,'Tasks');
 counter=0;

 %% retrieve output
 
 if verLessThan('matlab', '8')
    waitForState(job);
    out = getAllOutputArguments(job);
 else
    wait(job);
    out = fetchOutputs(job);
 end

 if options_.Cluster_settings == 3 % save files to localhost, if no shared filesystem is available
     for i=1:size(out,1)
         if (isfield(out{i},{'LocalFiles'}) )  
             for j=1:size(out{i}.LocalFiles,2)
                 if ((isfield(out{i}.LocalFiles(j),{'mat'}) ) == 1 ) && (isempty(cell2mat(out{i}.LocalFiles(j).mat)) == 0 )
                    fid=fopen([cell2mat(out{i}.LocalFiles(j).mat)],'w');
                    fwrite(fid,out{i}.LocalData(j).mat);
                    fclose(fid);
                 end
                 if ((isfield(out{i}.LocalFiles(j),{'log'}) ) == 1 ) && (isempty(cell2mat(out{i}.LocalFiles(j).log)) == 0 )              
                    fid=fopen([cell2mat(out{i}.LocalFiles(j).log)],'w');
                    fwrite(fid,out{i}.LocalData(j).log);
                    fclose(fid);
                 end
             end
         end
         out{i}.LocalFiles=[];
         out{i}.LocalData=[];
     end
 end
%% reformat output for dynare

for i=1:CPUsinUse
    fOutVar(i)=out{i};
end

% back to MasterParallel