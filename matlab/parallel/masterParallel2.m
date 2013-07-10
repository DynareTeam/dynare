function [fOutVar,nBlockPerCPU, CPUsinUse] = masterParallel2(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar,Parallel_info,initialize)

%dbstop in masterParallel2 at 118
%dbstop in random_walk_metropolis_hastings at 151

%% pack global variables, send to core function via input arguments

 globalVars = who('global');  
 
 for j=1:length(globalVars),
        eval(['global ',globalVars{j},';'])    
 end 
       
 for i=1:length(globalVars)
     name=globalVars{i};
     globalVars2.(globalVars{i})=eval(name);
 end
  
 fInputVar.global=globalVars2;
 fInputVar.Parallel = Parallel;
 
%a=random_walk_metropolis_hastings_core(fInputVar,1,2,0)
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
     
     case 0                                           % even distribution across CPUs
          BlocksPerCPU=NumberOfBlocks/totCPU;
          CPUsinUse=totCPU;
          
     case NumberOfBlocks                               % more CPUs than blocks, each block on a different CPU
          BlocksPerCPU=1;
          CPUsinUse=NumberOfBlocks;
     
     otherwise
          BlocksPerCPU=floor(NumberOfBlocks/totCPU);
          CPUsinUse=totCPU;
          remainder=rem(NumberOfBlocks,totCPU);        % the first |rem| workers get an extra-job
 end
 
 %% slice the task for workers
 slices=zeros(CPUsinUse,2);
 nBlockPerCPU=zeros(1,CPUsinUse);
 fblck=fBlock;

 for i=1:CPUsinUse
    nblck=(fblck+BlocksPerCPU)-1;
    if remainder && ( i <= remainder)                      % the first |rem| workers get an extra-job
       nblck=nblck+1;                 
    end
    slices(i,1)=fblck;
    slices(i,2)=nblck;
    nBlockPerCPU(1,i)=(nblck-fblck)+1;
    fblck=nblck+1;
 end
 
%% 
switch  options_.Cluster_settings
    
    case 2 % set up worker environment -- local
       
        parallel.defaultClusterProfile('local');
        cluster = parcluster();
        %cluster.ClusterMatlabRoot = '/opt/matlab2012b';                %Parallel(1).MatlabOctavePath;
        cluster.JobStorageLocation = [ '/jobdata/' Parallel(1).UserName];
        
        
    case 1 % set up cluster environment, maximum nuber of CPUs as defined in config-file
    
        if verLessThan('matlab', '8')
            cluster = findResource('scheduler', 'type', 'Torque');
            cluster.ClusterMatlabRoot = Parallel(1).MatlabOctavePath;
            cluster.ClusterSize = totCPU;
            cluster.ResourceTemplate = '-l nodes=^N^';
            cluster.DataLocation = [ '/jobdata/' Parallel(1).UserName];
        else
            cluster = parallel.cluster.Torque();
            
            cluster.ClusterMatlabRoot = '/opt/matlab2012b';                %Parallel(1).MatlabOctavePath;
            cluster.JobStorageLocation = [ '/jobdata/' Parallel(1).UserName];
            cluster.NumWorkers = totCPU;
        end
        
        cluster.HasSharedFilesystem = true;
        cluster.RcpCommand = 'scp';
        cluster.RshCommand = 'ssh';
        % You will get an email for EVERY SINGLE TASK.
        cluster.SubmitArguments = '-q verylong -l walltime=12:00:00 -m bea -M martin_alxander.westerberg@ecb.int';
        
    otherwise
        error('Unknown parallel option') 
end
 
 %% job-dispatch
 
 job = createJob(cluster);
 
  folders=regexp(path,':','split'); % current path settings of localhost
  folders(1,end+1)={cd};            % for dynamically created m-files
  set(job,'AdditionalPaths',folders)
  
 for i=1:CPUsinUse
    
     a=createTask(job, fname, 1, {fInputVar,slices(i,1),slices(i,2),0});  % WhoIam=0, we 'simulate' local computation
     a.CaptureDiary = true;
     
 end
 
 submit(job);
 get(job,'Tasks')
 
 %% retrieve output
 
 if verLessThan('matlab', '8')
    waitForState(job);
    out = getAllOutputArguments(job);
else
    wait(job);
    out = fetchOutputs(job);
 end
disp(out);

%% reformat output for dynare

for i=1:CPUsinUse
    fOutVar(i)=out{i};
end

% back to MasterParallel
