function [ ix2, ilogpo2, ModelName, MhDirectoryName, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = ...
    metropolis_hastings_initialization(TargetFun, xparam1, vv, mh_bounds,dataset_,options_,M_,estim_params_,bayestopt_,oo_)
%function [ ix2, ilogpo2, ModelName, MhDirectoryName, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = 
%    metropolis_hastings_initialization(TargetFun, xparam1, vv, mh_bounds, dataset_,options_,M_,estim_params_,bayestopt_,oo_)
% Metropolis-Hastings initialization.
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o dataset_              data structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%  
% OUTPUTS 
%   None  
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2011 Dynare Team
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

ix2 = [];
ilogpo2 = [];
ModelName = []; 
MhDirectoryName = [];
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
    ModelName = [M_.fname '_bvar'];
end

MhDirectoryName = CheckPath('metropolis',M_.dname);

nblck = options_.mh_nblck;
nruns = ones(nblck,1)*options_.mh_replic;
npar  = length(xparam1);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
d = chol(vv);

if ~options_.load_mh_file && ~options_.mh_recover
    %% Here we start a new metropolis-hastings, previous draws are not
    %% considered.
    if nblck > 1
        disp('MH: Multiple chains mode.')
    else
        disp('MH: One Chain mode.')
    end
    %% Delete old mh files...
    files = dir([ MhDirectoryName '/' ModelName '_mh*_blck*.mat']);
    if length(files)
        delete([ MhDirectoryName '/' ModelName '_mh*_blck*.mat']);
        disp('MH: Old _mh files successfully erased!')
    end
    file = dir([ MhDirectoryName '/metropolis.log']);
    if length(file)
        delete([ MhDirectoryName '/metropolis.log']);
        disp('MH: Old metropolis.log file successfully erased!')
        disp('MH: Creation of a new metropolis.log file.')
    end
    fidlog = fopen([MhDirectoryName '/metropolis.log'],'w');
    fprintf(fidlog,'%% MH log file (Dynare).\n');
    fprintf(fidlog,['%% ' datestr(now,0) '.\n']);
    fprintf(fidlog,' \n\n');
    fprintf(fidlog,'%% Session 1.\n');
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(nblck) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    %% Initial values...
    if nblck > 1% Case 1: multiple chains
        fprintf(fidlog,['  Initial values of the parameters:\n']);
        disp('MH: Searching for initial values...')
        ix2 = zeros(nblck,npar);
        ilogpo2 = zeros(nblck,1);
        for j=1:nblck
            validate    = 0;
            init_iter   = 0;
            trial = 1;
            while validate == 0 && trial <= 10
                candidate = rand_multivariate_normal( transpose(xparam1), d * options_.mh_init_scale, npar);
                if all(candidate(:) > mh_bounds(:,1)) && all(candidate(:) < mh_bounds(:,2)) 
                    ix2(j,:) = candidate;
                    ilogpo2(j) = - feval(TargetFun,ix2(j,:)',dataset_,options_,M_,estim_params_,bayestopt_,oo_);
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
                    disp(['MH: I couldn''t get a valid initial value in 100 trials.'])
                    disp(['MH: You should Reduce mh_init_scale...'])
                    disp(sprintf('MH: Parameter mh_init_scale is equal to %f.',options_.mh_init_scale))
                    options_.mh_init_scale = input('MH: Enter a new value...  ');
                    trial = trial+1;
                end
            end
            if trial > 10 && ~validate
                disp(['MH: I''m unable to find a starting value for block ' int2str(j)])
                return
            end
        end
        fprintf(fidlog,' \n');
        disp('MH: Initial values found!')
        disp(' ')
    else% Case 2: one chain (we start from the posterior mode)
        fprintf(fidlog,['  Initial values of the parameters:\n']);
        candidate = transpose(xparam1(:));%
        if all(candidate(:) > mh_bounds(:,1)) && all(candidate(:) < mh_bounds(:,2)) 
            ix2 = candidate;
            ilogpo2 = - feval(TargetFun,ix2',dataset_,options_,M_,estim_params_,bayestopt_,oo_);
            disp('MH: Initialization at the posterior mode.')
            disp(' ')
            fprintf(fidlog,['    Blck ' int2str(1) 'params:\n']);
            for i=1:length(ix2(1,:))
                fprintf(fidlog,['      ' int2str(i)  ':' num2str(ix2(1,i)) '\n']);
            end
            fprintf(fidlog,['    Blck ' int2str(1) 'logpo2:' num2str(ilogpo2) '\n']);
        else
            disp('MH: Initialization failed...')
            disp('MH: The posterior mode lies outside the prior bounds.')
            return
        end
        fprintf(fidlog,' \n');
    end
    fprintf(fidlog,' \n');
    fblck = 1;
    fline = ones(nblck,1);
    NewFile = ones(nblck,1);
    %%
    %% Creation of the mh-history file:
    %%
    file = dir([MhDirectoryName '/'  ModelName '_mh_history.mat']);
    if length(file)
        delete([ MhDirectoryName '/' ModelName '_mh_history.mat']);
        disp('MH: Old mh_history file successfully erased!')
    end
    AnticipatedNumberOfFiles = ceil(nruns(1)/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = nruns(1) - (AnticipatedNumberOfFiles-1)*MAX_nruns;
    record.Nblck = nblck;
    record.MhDraws = zeros(1,3);
    record.MhDraws(1,1) = nruns(1);
    record.MhDraws(1,2) = AnticipatedNumberOfFiles;
    record.MhDraws(1,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.AcceptationRates = zeros(1,nblck);
    % separate initializaton for each chain
    JSUM = 0;
    for j=1:nblck,
        % we set a different seed for the random generator for each block
        % then we record the corresponding random generator state (vector)
        set_dynare_seed(options_.DynareRandomStreams.seed+j);
        % record.Seeds keeps a vector of the random generator state and
        % not the scalar seed despite its name
        [record.Seeds(j).Unifor,record.Seeds(j).Normal] = get_dynare_random_generator_state();
    end
    record.InitialParameters = ix2;
    record.InitialLogLiK = ilogpo2;
    record.LastParameters = zeros(nblck,npar);
    record.LastLogPost = zeros(nblck,1);
    record.LastFileNumber = AnticipatedNumberOfFiles ;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    save([MhDirectoryName '/' ModelName '_mh_history.mat'],'record');  
    fprintf(fidlog,['  CREATION OF THE MH HISTORY FILE!\n\n']);
    fprintf(fidlog,['    Expected number of files per block.......: ' int2str(AnticipatedNumberOfFiles) '.\n']);
    fprintf(fidlog,['    Expected number of lines in the last file: ' ...
                    int2str(AnticipatedNumberOfLinesInTheLastFile) '.\n']);
    fprintf(fidlog,['\n']);
    for j = 1:nblck,
        fprintf(fidlog,['    Initial seed (randn) for chain number ',int2str(j),':\n']);
        for i=1:length(record.Seeds(j).Normal)
            fprintf(fidlog,['      ' num2str(record.Seeds(j).Normal(i)') '\n']);
        end
        fprintf(fidlog,['    Initial seed (rand) for chain number ',int2str(j),':\n']);
        for i=1:length(record.Seeds(j).Unifor)
            fprintf(fidlog,['      ' num2str(record.Seeds(j).Unifor(i)') '\n']);
        end
    end,
    fprintf(fidlog,' \n');
    fclose(fidlog);
elseif options_.load_mh_file && ~options_.mh_recover
    %% Here we consider previous mh files (previous mh did not crash).
    disp('MH: I''m loading past metropolis-hastings simulations...')
    file = dir([ MhDirectoryName '/'  ModelName '_mh_history.mat' ]);
    files = dir([ MhDirectoryName filesep ModelName '_mh*.mat']);
    if ~length(files)
        error('MH:: FAILURE! there is no MH file to load here!')
    end
    if ~length(file)
        error('MH:: FAILURE! there is no MH-history file!')
    else
        load([ MhDirectoryName '/'  ModelName '_mh_history.mat'])
    end
    fidlog = fopen([MhDirectoryName '/metropolis.log'],'a');
    fprintf(fidlog,'\n');
    fprintf(fidlog,['%% Session ' int2str(length(record.MhDraws(:,1))+1) '.\n']);
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(nblck) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    past_number_of_blocks = record.Nblck;
    if past_number_of_blocks ~= nblck
        disp('MH:: The specified number of blocks doesn''t match with the previous number of blocks!')
        disp(['MH:: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' ...
              int2str(past_number_of_blocks) '.'])
        disp(['MH:: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
        nblck = past_number_of_blocks;
        options_.mh_nblck = nblck;
    end
    % I read the last line of the last mh-file for initialization 
    % of the new metropolis-hastings simulations:
    LastFileNumber = record.LastFileNumber;
    LastLineNumber = record.LastLineNumber;
    if LastLineNumber < MAX_nruns
        NewFile = ones(nblck,1)*LastFileNumber;
    else
        NewFile = ones(nblck,1)*(LastFileNumber+1);
    end
    ilogpo2 = record.LastLogPost;
    ix2 = record.LastParameters;
    fblck = 1;
    fline = ones(nblck,1)*(LastLineNumber+1);
    NumberOfPreviousSimulations = sum(record.MhDraws(:,1),1);
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
    save([MhDirectoryName '/' ModelName '_mh_history.mat'],'record');
    disp(['MH: ... It''s done. I''ve loaded ' int2str(NumberOfPreviousSimulations) ' simulations.'])
    disp(' ')
    fclose(fidlog);
elseif options_.mh_recover
    %% The previous metropolis-hastings crashed before the end! I try to
    %% recover the saved draws...
    disp('MH: Recover mode!')
    disp(' ')
    file = dir([MhDirectoryName '/'  ModelName '_mh_history.mat']);
    if ~length(file)
        error('MH:: FAILURE! there is no MH-history file!')
    else
        load([ MhDirectoryName '/'  ModelName '_mh_history.mat'])
    end
    nblck = record.Nblck;% Number of "parallel" mcmc chains.
    options_.mh_nblck = nblck;
    if size(record.MhDraws,1) == 1
        OldMh = 0;% The crashed metropolis was the first session.
    else
        OldMh = 1;% The crashed metropolis wasn't the first session.
    end
    %% Default initialization:
    if OldMh
        ilogpo2 = record.LastLogPost;
        ix2 = record.LastParameters;
    else
        ilogpo2 = record.InitialLogLiK;
        ix2 = record.InitialParameters;
    end
    %% Set NewFile, a nblck*1 vector of integers, and fline (first line), a nblck*1 vector of integers.
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
    %% Set fblck (First block), an integer targeting the crashed mcmc chain.
    fblck = 1;
    %% How many mh files should we have ?
    ExpectedNumberOfMhFilesPerBlock = sum(record.MhDraws(:,2),1);
    ExpectedNumberOfMhFiles = ExpectedNumberOfMhFilesPerBlock*nblck;
    %% I count the total number of saved mh files...
    AllMhFiles = dir([MhDirectoryName '/' ModelName '_mh*_blck*.mat']);
    TotalNumberOfMhFiles = length(AllMhFiles);
    %% And I quit if I can't find a crashed mcmc chain. 
    if (TotalNumberOfMhFiles==ExpectedNumberOfMhFiles)
        disp('MH: It appears that you don''t need to use the mh_recover option!')
        disp('    You have to edit the mod file and remove the mh_recover') 
        disp('    in the estimation command')
        error
    end
    %% I count the number of saved mh files per block.
    NumberOfMhFilesPerBlock = zeros(nblck,1);
    for b = 1:nblck
        BlckMhFiles = dir([ MhDirectoryName '/' ModelName '_mh*_blck' int2str(b) '.mat']);
        NumberOfMhFilesPerBlock(b) = length(BlckMhFiles);
    end
    %% Is there a chain with less mh files than expected ? 
    while fblck <= nblck
        if  NumberOfMhFilesPerBlock(fblck) < ExpectedNumberOfMhFilesPerBlock
            disp(['MH: Chain ' int2str(fblck) ' is not complete!'])
            break
            % The mh_recover session will start from chain fblck.
        else
            disp(['MH: Chain ' int2str(fblck) ' is complete!'])
        end
        fblck = fblck+1;
    end
    %% How many mh-files are saved in this block?
    NumberOfSavedMhFilesInTheCrashedBlck = NumberOfMhFilesPerBlock(fblck);
    %% Correct the number of saved mh files if the crashed metropolis was not the first session (so
    %% that NumberOfSavedMhFilesInTheCrashedBlck is the number of saved mh files in the crashed chain 
    %% of the current session).  
    if OldMh
        NumberOfSavedMhFilesInTheCrashedBlck = NumberOfSavedMhFilesInTheCrashedBlck - LastFileNumberInThePreviousMh;
    end
    NumberOfSavedMhFiles = NumberOfSavedMhFilesInTheCrashedBlck+LastFileNumberInThePreviousMh;
    %% Correct initial conditions.
    if NumberOfSavedMhFiles
        load([MhDirectoryName '/' ModelName '_mh' int2str(NumberOfSavedMhFiles) '_blck' int2str(fblck) '.mat']);
        ilogpo2(fblck) = logpo2(end);
        ix2(fblck,:) = x2(end,:);
    end
end