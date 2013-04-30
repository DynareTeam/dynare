function myoutput = random_walk_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% This function contain the most computationally intensive portion of code in
% random_walk_metropolis_hastings (the 'for xxx = fblck:nblck' loop). The branches in 'for'
% cycle and are completely independent than suitable to be executed in parallel way.
%
% INPUTS
%   o myimput            [struc]     The mandatory variables for local/remote
%                                    parallel computing obtained from random_walk_metropolis_hastings.m
%                                    function.
%   o fblck and nblck    [integer]   The Metropolis-Hastings chains.
%   o whoiam             [integer]   In concurrent programming a modality to refer to the differents thread running in parallel is needed.
%                                    The integer whoaim is the integer that
%                                    allows us to distinguish between them. Then it is the index number of this CPU among all CPUs in the
%                                    cluster.
%   o ThisMatlab         [integer]   Allows us to distinguish between the
%                                    'main' matlab, the slave matlab worker, local matlab, remote matlab,
%                                     ... Then it is the index number of this slave machine in the cluster.
% OUTPUTS
%   o myoutput  [struc]
%               If executed without parallel is the original output of 'for b =
%               fblck:nblck' otherwise a portion of it computed on a specific core or
%               remote machine. In this case:
%                               record;
%                               irun;
%                               NewFile;
%                               OutputFileName
%
% ALGORITHM
%   Portion of Metropolis-Hastings.
%
% SPECIAL REQUIREMENTS.
%   None.

% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code sutable to be executed in parallel on multi core or cluster machine,
% is removed from this function and placed in random_walk_metropolis_hastings_core.m funtion.
% Then the DYNARE parallel package contain a set of pairs matlab functios that can be executed in
% parallel and called name_function.m and name_function_core.m.
% In addition in the parallel package we have second set of functions used
% to manage the parallel computation.
%
% This function was the first function to be parallelized, later other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management funtions.


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

if nargin<4,
    whoiam=0;
end

% reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

TargetFun=myinputs.TargetFun;
ProposalFun=myinputs.ProposalFun;
xparam1=myinputs.xparam1;
vv=myinputs.vv;
mh_bounds=myinputs.mh_bounds;
ix2=myinputs.ix2;
ilogpo2=myinputs.ilogpo2;
ModelName=myinputs.ModelName;
fline=myinputs.fline;
npar=myinputs.npar;
nruns=myinputs.nruns;
NewFile=myinputs.NewFile;
MAX_nruns=myinputs.MAX_nruns;
d=myinputs.d;
InitSizeArray=myinputs.InitSizeArray;
record=myinputs.record;
dataset_ = myinputs.dataset_;
bayestopt_ = myinputs.bayestopt_;
estim_params_ = myinputs.estim_params_;
options_ = myinputs.options_;
M_ = myinputs.M_;
oo_ = myinputs.oo_;
varargin=myinputs.varargin;

% Necessary only for remote computing!
if whoiam
    Parallel=myinputs.Parallel;
    % initialize persistent variables in priordens()
    priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, ...
        bayestopt_.p3,bayestopt_.p4,1);
end

MhDirectoryName = CheckPath('metropolis',M_.dname);

options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);
if strcmpi(ProposalFun,'rand_multivariate_normal')
    n = npar;
elseif strcmpi(ProposalFun,'rand_multivariate_student')
    n = options_.student_degrees_of_freedom;
end

%%%%
%%%% NOW i run the (nblck-fblck+1) metropolis-hastings chains
%%%%

proposal_covariance_Cholesky_decomposition = d*diag(bayestopt_.jscale);

jloop=0;

JSUM = 0;
for b = fblck:nblck,
    jloop=jloop+1;
    try
        % this will not work if the master uses a random generator not
        % available in the slave (different Matlab version or
        % Matlab/Octave cluster). Therefor the trap.
        
        % this set the random generator type (the seed is useless but
        % needed by the function)
        set_dynare_seed(options_.DynareRandomStreams.algo,...
                        options_.DynareRandomStreams.seed);
        % this set the state 
        set_dynare_random_generator_state(record.Seeds(b).Unifor, ...
                                          record.Seeds(b).Normal);
    catch
        % if the state set by master is incompatible with the slave, we
        % only reseed 
        set_dynare_seed(options_.DynareRandomStreams.seed+b);
    end
    if (options_.load_mh_file~=0) && (fline(b)>1) && OpenOldFile(b)
        load([pwd filesep MhDirectoryName filesep ModelName '_mh' int2str(NewFile(b)) ...
            '_blck' int2str(b) '.mat'])
        x2 = [x2;zeros(InitSizeArray(b)-fline(b)+1,npar)];
        logpo2 = [logpo2;zeros(InitSizeArray(b)-fline(b)+1,1)];
        OpenOldFile(b) = 0;
    else
        x2 = zeros(InitSizeArray(b),npar);
        logpo2 = zeros(InitSizeArray(b),1);
    end
    if whoiam
        prc0=(b-fblck)/(nblck-fblck+1)*(exist('OCTAVE_VERSION') || options_.console_mode);
        hh = dyn_waitbar({prc0,whoiam,options_.parallel(ThisMatlab)},['MH (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
    else
        hh = dyn_waitbar(0,['Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
        set(hh,'Name','Metropolis-Hastings');
    end
    isux = 0;
    jsux = 0;
    irun = fline(b);
    j = 1;
    while j <= nruns(b)
        par = feval(ProposalFun, ix2(b,:), proposal_covariance_Cholesky_decomposition, n);
        if all( par(:) > mh_bounds(:,1) ) && all( par(:) < mh_bounds(:,2) )
            try
                logpost = - feval(TargetFun, par(:),dataset_,options_,M_,estim_params_,bayestopt_,oo_);
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
        if (logpost > -inf) && (log(rand) < logpost-ilogpo2(b))
            x2(irun,:) = par;
            ix2(b,:) = par;
            logpo2(irun) = logpost;
            ilogpo2(b) = logpost;
            isux = isux + 1;
            jsux = jsux + 1;
        else
            x2(irun,:) = ix2(b,:);
            logpo2(irun) = ilogpo2(b);
        end
        prtfrc = j/nruns(b);
%         if exist('OCTAVE_VERSION') || options_.console_mode
%             if mod(j, 10) == 0
%                 if exist('OCTAVE_VERSION')
%                     if (whoiam==0)
%                         printf('MH: Computing Metropolis-Hastings (chain %d/%d): %3.f%% done, acception rate: %3.f%%\r', b, nblck, 100 * prtfrc, 100 * isux / j);
%                     end
%                 else
%                     s0=repmat('\b',1,length(newString));
%                     newString=sprintf('MH: Computing Metropolis-Hastings (chain %d/%d): %3.f%% done, acceptance rate: %3.f%%', b, nblck, 100 * prtfrc, 100 * isux / j);
%                     fprintf([s0,'%s'],newString);
%                 end
%             end
%             if mod(j,50)==0 && whoiam
%                 %             keyboard;
%                 if (strcmp([options_.parallel(ThisMatlab).MatlabOctavePath], 'octave'))
%                     waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) '), ' sprintf('accept. %3.f%%',100 * isux / j)];
%                     fMessageStatus(prtfrc,whoiam,waitbarString, waitbarTitle, options_.parallel(ThisMatlab));
%                 else
%                     waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) '), ' sprintf('accept. %3.f%%', 100 * isux/j)];
%                     fMessageStatus((b-fblck)/(nblck-fblck+1)+prtfrc/(nblck-fblck+1),whoiam,waitbarString, '', options_.parallel(ThisMatlab));
%                 end
%             end
%         else
%             if mod(j, 3)==0 && ~whoiam
%                 waitbar(prtfrc,hh,[ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)]);
%             elseif mod(j,50)==0 && whoiam,
%                 %             keyboard;
%                 waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)];
%                 fMessageStatus(prtfrc,whoiam,waitbarString, waitbarTitle, options_.parallel(ThisMatlab));
%             end
%         end
        if (mod(j, 3)==0 && ~whoiam) || (mod(j,50)==0 && whoiam)
            dyn_waitbar(prtfrc,hh,[ 'MH (' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('acceptation rate %4.3f', isux/j)]);
        end
        if (irun == InitSizeArray(b)) || (j == nruns(b)) % Now I save the simulations
            save([MhDirectoryName '/' ModelName '_mh' int2str(NewFile(b)) '_blck' int2str(b) '.mat'],'x2','logpo2');
            fidlog = fopen([MhDirectoryName '/metropolis.log'],'a');
            fprintf(fidlog,['\n']);
            fprintf(fidlog,['%% Mh' int2str(NewFile(b)) 'Blck' int2str(b) ' (' datestr(now,0) ')\n']);
            fprintf(fidlog,' \n');
            fprintf(fidlog,['  Number of simulations.: ' int2str(length(logpo2)) '\n']);
            fprintf(fidlog,['  Acceptation rate......: ' num2str(jsux/length(logpo2)) '\n']);
            fprintf(fidlog,['  Posterior mean........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(mean(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(mean(logpo2)) '\n']);
            fprintf(fidlog,['  Minimum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(min(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(min(logpo2)) '\n']);
            fprintf(fidlog,['  Maximum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(max(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(max(logpo2)) '\n']);
            fprintf(fidlog,' \n');
            fclose(fidlog);
            jsux = 0;
            if j == nruns(b) % I record the last draw...
                record.LastParameters(b,:) = x2(end,:);
                record.LastLogPost(b) = logpo2(end);
            end
            % size of next file in chain b
            InitSizeArray(b) = min(nruns(b)-j,MAX_nruns);
            % initialization of next file if necessary
            if InitSizeArray(b)
                x2 = zeros(InitSizeArray(b),npar);
                logpo2 = zeros(InitSizeArray(b),1);
                NewFile(b) = NewFile(b) + 1;
                irun = 0;
            end
        end
        j=j+1;
        irun = irun + 1;
    end% End of the simulations for one mh-block.
    record.AcceptationRates(b) = isux/j;
%     if exist('OCTAVE_VERSION') || options_.console_mode || whoiam
%         if exist('OCTAVE_VERSION')
%             printf('\n');
%         else
%             fprintf('\n');
%         end
%         diary on;
%     else %if ~whoiam
%         close(hh);
%     end
    dyn_waitbar_close(hh);
    [record.Seeds(b).Unifor, record.Seeds(b).Normal] = get_dynare_random_generator_state();
    OutputFileName(jloop,:) = {[MhDirectoryName,filesep], [ModelName '_mh*_blck' int2str(b) '.mat']};
end% End of the loop over the mh-blocks.


myoutput.record = record;
myoutput.irun = irun;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;