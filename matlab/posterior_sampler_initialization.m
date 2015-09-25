function [ ix2, ilogpo2, ModelName, MetropolisFolder, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = ...
    posterior_sampler_initialization(TargetFun, xparam1, vv, mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
%function [ ix2, ilogpo2, ModelName, MetropolisFolder, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = ...
%     metropolis_hastings_initialization(TargetFun, xparam1, vv, mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% Metropolis-Hastings initialization.
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o dataset_              data structure
%   o dataset_info          dataset info structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%  
% OUTPUTS 
%   o ix2                   [double]   (nblck*npar) vector of starting points for different chains
%   o ilogpo2               [double]   (nblck*1) vector of initial posterior values for different chains
%   o ModelName             [string]    name of the mod-file
%   o MetropolisFolder      [string]    path to the Metropolis subfolder
%   o fblck                 [scalar]    number of the first MH chain to be run (not equal to 1 in case of recovery)   
%   o fline                 [double]   (nblck*1) vector of first draw in each chain (not equal to 1 in case of recovery)   
%   o npar                  [scalar]    number of parameters estimated
%   o nblck                 [scalar]    Number of MCM chains requested   
%   o nruns                 [double]   (nblck*1) number of draws in each chain 
%   o NewFile               [scalar]    (nblck*1) vector storing the number
%                                       of the first MH-file to created for each chain when saving additional
%                                       draws
%   o MAX_nruns             [scalar]    maximum number of draws in each MH-file on the harddisk
%   o d                     [double]   (p*p) matrix, Cholesky decomposition of the posterior covariance matrix (at the mode).
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2015 Dynare Team
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

%Initialize outputs
ix2 = [];
ilogpo2 = [];
ModelName = []; 
MetropolisFolder = [];
fblck = [];
fline = [];
npar = [];
nblck = [];
nruns = [];
NewFile = [];
MAX_nruns = [];
d = [];

ModelName = M_.fname;
if ~isempty(M_.bvar)
    ModelName = [ModelName '_bvar'];
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
BaseName = [MetropolisFolder filesep ModelName];

nblck = options_.mh_nblck;
nruns = ones(nblck,1)*options_.mh_replic;
npar  = length(xparam1);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
d = chol(vv);

if ~options_.load_mh_file && ~options_.mh_recover
    % Here we start a new Metropolis-Hastings, previous draws are discarded.
    if nblck > 1
        disp('Estimation::mcmc: Multiple chains mode.')
    else
        disp('Estimation::mcmc: One Chain mode.')
    end
    % Delete old mh files if any...
    files = dir([BaseName '_mh*_blck*.mat']);
    if length(files)
        delete([BaseName '_mh*_blck*.mat']);
        disp('Estimation::mcmc: Old mh-files successfully erased!')
    end
    % Delete old Metropolis log file.
    file = dir([ MetropolisFolder '/metropolis.log']);
    if length(file)
        delete([ MetropolisFolder '/metropolis.log']);
        disp('Estimation::mcmc: Old metropolis.log file successfully erased!')
        disp('Estimation::mcmc: Creation of a new metropolis.log file.')
    end
    fidlog = fopen([MetropolisFolder '/metropolis.log'],'w');
    fprintf(fidlog,'%% MH log file (Dynare).\n');
    fprintf(fidlog,['%% ' datestr(now,0) '.\n']);
    fprintf(fidlog,' \n\n');
    fprintf(fidlog,'%% Session 1.\n');
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(nblck) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    % Find initial values for the nblck chains...
    if isempty(d),
        prior_draw(1);
    end
    if nblck > 1% Case 1: multiple chains
        fprintf(fidlog,['  Initial values of the parameters:\n']);
        disp('Estimation::mcmc: Searching for initial values...')
        ix2 = zeros(nblck,npar);
        ilogpo2 = zeros(nblck,1);
        for j=1:nblck
            validate    = 0;
            init_iter   = 0;
            trial = 1;
            while validate == 0 && trial <= 10
                if isempty(d),
                    candidate = prior_draw(0);
                else
                    candidate = rand_multivariate_normal( transpose(xparam1), d * options_.mh_init_scale, npar);
                end
                if all(candidate(:) > mh_bounds.lb) && all(candidate(:) < mh_bounds.ub) 
                    ix2(j,:) = candidate;
                    ilogpo2(j) = - feval(TargetFun,ix2(j,:)',dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
                    if ~isfinite(ilogpo2(j)) % if returned log-density is
                                             % Inf or Nan (penalized value)
                        validate = 0;
                    else
                        fprintf(fidlog,['    Blck ' int2str(j) ':\n']);
                        for i=1:length(ix2(1,:))
                            fprintf(fidlog,['      params:' int2str(i) ': ' num2str(ix2(j,i)) '\n']);
                        end
                        fprintf(fidlog,['      logpo2: ' num2str(ilogpo2(j)) '\n']);
                        j = j+1;
                        validate = 1;
                    end
                end
                init_iter = init_iter + 1;
                if init_iter > 100 && validate == 0
                    disp(['Estimation::mcmc: I couldn''t get a valid initial value in 100 trials.'])
                    if options_.nointeractive
                        disp(['Estimation::mcmc: I reduce mh_init_scale by 10 percent:'])
                        options_.mh_init_scale = .9*options_.mh_init_scale;
                        disp(sprintf('Estimation::mcmc: Parameter mh_init_scale is now equal to %f.',options_.mh_init_scale))
                    else
                        disp(['Estimation::mcmc: You should reduce mh_init_scale...'])
                        disp(sprintf('Estimation::mcmc: Parameter mh_init_scale is equal to %f.',options_.mh_init_scale))
                        options_.mh_init_scale = input('Estimation::mcmc: Enter a new value...  ');
                    end
                    trial = trial+1;
                end
            end
            if trial > 10 && ~validate
                disp(['Estimation::mcmc: I''m unable to find a starting value for block ' int2str(j)])
                return
            end
        end
        fprintf(fidlog,' \n');
        disp('Estimation::mcmc: Initial values found!')
        skipline()
    else% Case 2: one chain (we start from the posterior mode)
        fprintf(fidlog,['  Initial values of the parameters:\n']);
        candidate = transpose(xparam1(:));%
        if all(candidate(:) > mh_bounds.lb) && all(candidate(:) < mh_bounds.ub) 
            ix2 = candidate;
            ilogpo2 = - feval(TargetFun,ix2',dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
            disp('Estimation::mcmc: Initialization at the posterior mode.')
            skipline()
            fprintf(fidlog,['    Blck ' int2str(1) 'params:\n']);
            for i=1:length(ix2(1,:))
                fprintf(fidlog,['      ' int2str(i)  ':' num2str(ix2(1,i)) '\n']);
            end
            fprintf(fidlog,['    Blck ' int2str(1) 'logpo2:' num2str(ilogpo2) '\n']);
        else
            disp('Estimation::mcmc: Initialization failed...')
            disp('Estimation::mcmc: The posterior mode lies outside the prior bounds.')
            return
        end
        fprintf(fidlog,' \n');
    end
    fprintf(fidlog,' \n');
    fblck = 1;
    fline = ones(nblck,1);
    NewFile = ones(nblck,1);
    % Delete the mh-history files
    delete_mh_history_files(MetropolisFolder, ModelName);
    %  Create a new record structure
    fprintf(['Estimation::mcmc: Write details about the MCMC... ']);
    AnticipatedNumberOfFiles = ceil(nruns(1)/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = nruns(1) - (AnticipatedNumberOfFiles-1)*MAX_nruns;
    record.Sampler = options_.posterior_sampling_method;
    record.Nblck = nblck;
    record.MhDraws = zeros(1,3);
    record.MhDraws(1,1) = nruns(1);
    record.MhDraws(1,2) = AnticipatedNumberOfFiles;
    record.MhDraws(1,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.AcceptanceRatio = zeros(1,nblck);
    record.FunctionEvalPerIteration = NaN(1,nblck);
    for j=1:nblck
        % we set a different seed for the random generator for each block then we record the corresponding random generator state (vector)
        set_dynare_seed(options_.DynareRandomStreams.seed+j);
        % record.Seeds keeps a vector of the random generator state and not the scalar seed despite its name
        [record.InitialSeeds(j).Unifor,record.InitialSeeds(j).Normal] = get_dynare_random_generator_state();
    end
    record.InitialParameters = ix2;
    record.InitialLogPost = ilogpo2;
    record.LastParameters = zeros(nblck,npar);
    record.LastLogPost = zeros(nblck,1);
    record.LastFileNumber = AnticipatedNumberOfFiles ;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    record.MCMCConcludedSuccessfully = 0;
    fprintf('Ok!\n');
    id = write_mh_history_file(MetropolisFolder, ModelName, record);
    disp(['Estimation::mcmc: Details about the MCMC are available in ' BaseName '_mh_history_' num2str(id) '.mat'])
    skipline()
    fprintf(fidlog,['  CREATION OF THE MH HISTORY FILE!\n\n']);
    fprintf(fidlog,['    Expected number of files per block.......: ' int2str(AnticipatedNumberOfFiles) '.\n']);
    fprintf(fidlog,['    Expected number of lines in the last file: ' int2str(AnticipatedNumberOfLinesInTheLastFile) '.\n']);
    fprintf(fidlog,['\n']);
    for j = 1:nblck,
        fprintf(fidlog,['    Initial state of the Gaussian random number generator for chain number ',int2str(j),':\n']);
        for i=1:length(record.InitialSeeds(j).Normal)
            fprintf(fidlog,['      ' num2str(record.InitialSeeds(j).Normal(i)') '\n']);
        end
        fprintf(fidlog,['    Initial state of the Uniform random number generator for chain number ',int2str(j),':\n']);
        for i=1:length(record.InitialSeeds(j).Unifor)
            fprintf(fidlog,['      ' num2str(record.InitialSeeds(j).Unifor(i)') '\n']);
        end
    end,
    fprintf(fidlog,' \n');
    fclose(fidlog);
elseif options_.load_mh_file && ~options_.mh_recover
    % Here we consider previous mh files (previous mh did not crash).
    disp('Estimation::mcmc: I am loading past Metropolis-Hastings simulations...')
    load_last_mh_history_file(MetropolisFolder, ModelName);
    if ~isnan(record.MCMCConcludedSuccessfully) && ~record.MCMCConcludedSuccessfully
        error('Estimation::mcmc: You are trying to load an MCMC that did not finish successfully. Please use mh_recover.')
    end
    record.MCMCConcludedSuccessfully=0; %reset indicator for this run
    mh_files = dir([ MetropolisFolder filesep ModelName '_mh*.mat']);
    if ~length(mh_files)
        error('Estimation::mcmc: I cannot find any MH file to load here!')
    end
    fidlog = fopen([MetropolisFolder '/metropolis.log'],'a');
    fprintf(fidlog,'\n');
    fprintf(fidlog,['%% Session ' int2str(length(record.MhDraws(:,1))+1) '.\n']);
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(nblck) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    past_number_of_blocks = record.Nblck;
    if past_number_of_blocks ~= nblck
        disp('Estimation::mcmc: The specified number of blocks doesn''t match with the previous number of blocks!')
        disp(['Estimation::mcmc: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' int2str(past_number_of_blocks) '.'])
        disp(['Estimation::mcmc: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
        nblck = past_number_of_blocks;
        options_.mh_nblck = nblck;
    end
    % I read the last line of the last mh-file for initialization of the new Metropolis-Hastings simulations:
    LastFileNumber = record.LastFileNumber;
    LastLineNumber = record.LastLineNumber;
    if LastLineNumber < MAX_nruns
        NewFile = ones(nblck,1)*LastFileNumber;
        fline = ones(nblck,1)*(LastLineNumber+1);
    else
        NewFile = ones(nblck,1)*(LastFileNumber+1);
        fline = ones(nblck,1);
    end
    ilogpo2 = record.LastLogPost;
    ix2 = record.LastParameters;
    fblck = 1;
    NumberOfPreviousSimulations = sum(record.MhDraws(:,1),1);
    fprintf('Estimation::mcmc: I am writing a new mh-history file... ');
    record.MhDraws = [record.MhDraws;zeros(1,3)];
    NumberOfDrawsWrittenInThePastLastFile = MAX_nruns - LastLineNumber;
    NumberOfDrawsToBeSaved = nruns(1) - NumberOfDrawsWrittenInThePastLastFile;
    AnticipatedNumberOfFiles = ceil(NumberOfDrawsToBeSaved/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = NumberOfDrawsToBeSaved - (AnticipatedNumberOfFiles-1)*MAX_nruns;  
    record.LastFileNumber = LastFileNumber + AnticipatedNumberOfFiles;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    record.MhDraws(end,1) = nruns(1);
    record.MhDraws(end,2) = AnticipatedNumberOfFiles;
    record.MhDraws(end,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.InitialSeeds = record.LastSeeds;
    write_mh_history_file(MetropolisFolder, ModelName, record);
    fprintf('Done.\n')
    disp(['Estimation::mcmc: Ok. I have loaded ' int2str(NumberOfPreviousSimulations) ' simulations.'])
    skipline()
    fclose(fidlog);
elseif options_.mh_recover
    % The previous metropolis-hastings crashed before the end! I try to recover the saved draws...
    disp('Estimation::mcmc: Recover mode!')
    load_last_mh_history_file(MetropolisFolder, ModelName);
    nblck = record.Nblck;% Number of "parallel" mcmc chains.
    options_.mh_nblck = nblck;
    if size(record.MhDraws,1) == 1
        OldMh = 0;% The crashed metropolis was the first session.
    else
        OldMh = 1;% The crashed metropolis wasn't the first session.
    end
    % Default initialization:
    if OldMh
        ilogpo2 = record.LastLogPost;
        ix2 = record.LastParameters;
    else
        ilogpo2 = record.InitialLogPost;
        ix2 = record.InitialParameters;
    end
    % Set NewFile, a nblck*1 vector of integers, and fline (first line), a nblck*1 vector of integers.
    if OldMh
        LastLineNumberInThePreviousMh = record.MhDraws(end-1,3);% Number of lines in the last mh files of the previous session.
        LastFileNumberInThePreviousMh = sum(record.MhDraws(1:end-1,2),1);% Number of mh files in the the previous sessions.
        if LastLineNumberInThePreviousMh < MAX_nruns% Test if the last mh files of the previous session are not complete (yes)
            NewFile = ones(nblck,1)*LastFileNumberInThePreviousMh;
            fline = ones(nblck,1)*(LastLineNumberInThePreviousMh+1);
        else% The last mh files of the previous session are complete.
            NewFile = ones(nblck,1)*(LastFileNumberInThePreviousMh+1);
            fline = ones(nblck,1);
        end
    else
        LastLineNumberInThePreviousMh = 0;
        LastFileNumberInThePreviousMh = 0;
        NewFile = ones(nblck,1);
        fline = ones(nblck,1);
    end
    % Set fblck (First block), an integer targeting the crashed mcmc chain.
    fblck = 1;
    % How many mh files should we have ?
    ExpectedNumberOfMhFilesPerBlock = sum(record.MhDraws(:,2),1);
    ExpectedNumberOfMhFiles = ExpectedNumberOfMhFilesPerBlock*nblck;
    % I count the total number of saved mh files...
    AllMhFiles = dir([BaseName '_mh*_blck*.mat']);
    TotalNumberOfMhFiles = length(AllMhFiles);
    % And I quit if I can't find a crashed mcmc chain. 
    if (TotalNumberOfMhFiles==ExpectedNumberOfMhFiles)
        disp('Estimation::mcmc: It appears that you don''t need to use the mh_recover option!')
        disp('                  You have to edit the mod file and remove the mh_recover option') 
        disp('                  in the estimation command')
        error('Estimation::mcmc: mh_recover option not required!')
    end
    % I count the number of saved mh files per block.
    NumberOfMhFilesPerBlock = zeros(nblck,1);
    for b = 1:nblck
        BlckMhFiles = dir([BaseName '_mh*_blck' int2str(b) '.mat']);
        NumberOfMhFilesPerBlock(b) = length(BlckMhFiles);
    end
    % Is there a chain with less mh files than expected ? 
    while fblck <= nblck
        if  NumberOfMhFilesPerBlock(fblck) < ExpectedNumberOfMhFilesPerBlock
            disp(['Estimation::mcmc: Chain ' int2str(fblck) ' is not complete!'])
            break
            % The mh_recover session will start from chain fblck.
        else
            disp(['Estimation::mcmc: Chain ' int2str(fblck) ' is complete!'])
        end
        fblck = fblck+1;
    end
    % How many mh-files are saved in this block?
    NumberOfSavedMhFilesInTheCrashedBlck = NumberOfMhFilesPerBlock(fblck);
    % Correct the number of saved mh files if the crashed Metropolis was not the first session (so
    % that NumberOfSavedMhFilesInTheCrashedBlck is the number of saved mh files in the crashed chain 
    % of the current session).  
    if OldMh
        NumberOfSavedMhFilesInTheCrashedBlck = NumberOfSavedMhFilesInTheCrashedBlck - LastFileNumberInThePreviousMh;
    end
    NumberOfSavedMhFiles = NumberOfSavedMhFilesInTheCrashedBlck+LastFileNumberInThePreviousMh;
    % Correct initial conditions.
    if NumberOfSavedMhFiles
        load([BaseName '_mh' int2str(NumberOfSavedMhFiles) '_blck' int2str(fblck) '.mat']);
        ilogpo2(fblck) = logpo2(end);
        ix2(fblck,:) = x2(end,:);
    end
end