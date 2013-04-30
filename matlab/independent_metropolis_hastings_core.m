function myoutput = independent_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% The most computationally intensive portion of code in
% independent_metropolis_hastings (the 'for xxx = fblck:nblck' cycle).
% See the comment in random_walk_metropolis_hastings_core.m funtion.
%
% INPUTS
%   See See the comment in random_walk_metropolis_hastings_core.m funtion.

% OUTPUTS
%   See See the comment in random_walk_metropolis_hastings_core.m funtion.
%
% ALGORITHM
%   Portion of Independing Metropolis-Hastings.
%
% SPECIAL REQUIREMENTS.
%   None.
%
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


global bayestopt_ estim_params_ options_  M_ oo_ objective_function_penalty_base

% Reshape 'myinputs' for local computation.
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
varargin=myinputs.varargin;

if whoiam
    Parallel=myinputs.Parallel;
    % initialize persistent variables in priordens()
    priordens(xparam1',bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, ...
              bayestopt_.p3,bayestopt_.p4,1);
end

% (re)Set the penalty.
objective_function_penalty_base = Inf;

MhDirectoryName = CheckPath('metropolis',M_.dname);

OpenOldFile = ones(nblck,1);
if strcmpi(ProposalFun,'rand_multivariate_normal')
    n = npar;
    ProposalDensity = 'multivariate_normal_pdf';
elseif strcmpi(ProposalFun,'rand_multivariate_student')
    n = options_.student_degrees_of_freedom;
    ProposalDensity = 'multivariate_student_pdf';
end
% load([MhDirectoryName '/' ModelName '_mh_history.mat'],'record');
%%%%
%%%% NOW i run the (nblck-fblck+1) metropolis-hastings chains
%%%%

if any(isnan(bayestopt_.jscale))
    if exist([ModelName '_optimal_mh_scale_parameter.mat'])% This file is created by mode_compute=6.
        load([ModelName '_optimal_mh_scale_parameter'])
        proposal_covariance = d*Scale;
    else
        error('mh:: Something is wrong. I can''t figure out the value of the scale parameter.')
    end
else
    proposal_covariance = d*diag(bayestopt_.jscale);
end

jloop=0;

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
    if (options_.load_mh_file~=0)  && (fline(b)>1) && OpenOldFile(b)
        load([pwd filesep MhDirectoryName filesep ModelName '_mh' int2str(NewFile(b)) ...
              '_blck' int2str(b) '.mat'])
        x2 = [x2;zeros(InitSizeArray(b)-fline(b)+1,npar)];
        logpo2 = [logpo2;zeros(InitSizeArray(b)-fline(b)+1,1)];
        OpenOldFile(b) = 0;
    else
        x2 = zeros(InitSizeArray(b),npar);
        logpo2 = zeros(InitSizeArray(b),1);
    end
    if exist('OCTAVE_VERSION') || options_.console_mode
        diary off
        disp(' ')
    elseif whoiam
        %       keyboard;
        waitbarString = ['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...'];
        %       waitbarTitle=['Metropolis-Hastings ',options_.parallel(ThisMatlab).ComputerName];
        if options_.parallel(ThisMatlab).Local,
            waitbarTitle=['Local '];
        else
            waitbarTitle=[options_.parallel(ThisMatlab).ComputerName];
        end
        fMessageStatus(0,whoiam,waitbarString, waitbarTitle, options_.parallel(ThisMatlab));
    else,
        hh = waitbar(0,['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
        set(hh,'Name','Metropolis-Hastings');

    end
    isux = 0;
    jsux = 0;
    irun = fline(b);
    j = 1;
    while j <= nruns(b)
        par = feval(ProposalFun, xparam1, proposal_covariance, n);
        if all( par(:) > mh_bounds(:,1) ) && all( par(:) < mh_bounds(:,2) )
            try
                logpost = - feval(TargetFun, par(:),varargin{:});
            catch,
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
        r = logpost - ilogpo2(b) + ...
            log(feval(ProposalDensity, ix2(b,:), xparam1, proposal_covariance, n)) - ...
            log(feval(ProposalDensity, par, xparam1, proposal_covariance, n));
        if (logpost > -inf) && (log(rand) < r)
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
        if exist('OCTAVE_VERSION') || options_.console_mode
            if mod(j, 10) == 0
                if exist('OCTAVE_VERSION')
                    if (whoiam==0),
                        printf('MH: Computing Metropolis-Hastings (chain %d/%d): %3.f%% done, acception rate: %3.f%%\r', b, nblck, 100 * prtfrc, 100 * isux / j);
                    end
                else
                    fprintf('   MH: Computing Metropolis-Hastings (chain %d/%d): %3.f \b%% done, acception rate: %3.f \b%%\r', b, nblck, 100 * prtfrc, 100 * isux / j);
                end
            end
            if mod(j,50)==0 && whoiam,
                %             keyboard;
                waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) '), ' sprintf('accept. %3.f%%%%', 100 * isux/j)];
                fMessageStatus(prtfrc,whoiam,waitbarString, '', options_.parallel(ThisMatlab))
            end
        else
            if mod(j, 3)==0 && ~whoiam
                waitbar(prtfrc,hh,[ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)]);
            elseif mod(j,50)==0 && whoiam,
                %             keyboard;
                waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)];
                fMessageStatus(prtfrc,whoiam,waitbarString, waitbarTitle, options_.parallel(ThisMatlab))
            end
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
            fprintf(fidlog,['  Minimum value.........:\n']);;
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
    if exist('OCTAVE_VERSION') || options_.console_mode
        printf('\n');
        diary on
    elseif ~whoiam
        close(hh);
    end
    [record.Seeds(b).Unifor, record.Seeds(b).Normal] = get_dynare_random_generator_state();
    OutputFileName(jloop,:) = {[MhDirectoryName,filesep], [ModelName '_mh*_blck' int2str(b) '.mat']};
end% End of the loop over the mh-blocks.


myoutput.record = record;
myoutput.irun = irun;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;


