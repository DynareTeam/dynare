function prior_posterior_statistics(type,dataset)

% function prior_posterior_statistics(type,dataset)
% Computes Monte Carlo filter smoother and forecasts
%
% INPUTS
%    type:         posterior
%                  prior
%                  gsa
%    dataset:      data structure
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% PARALLEL CONTEXT
% See the comments random_walk_metropolis_hastings.m funtion.


% Copyright (C) 2005-2013 Dynare Team
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

global options_ estim_params_ oo_ M_ bayestopt_

localVars=[];

Y = dataset.data;
gend = dataset.info.ntobs;
data_index = dataset.missing.aindex;
missing_value = dataset.missing.state;
bayestopt_.mean_varobs = dataset.descriptive.mean';

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
naK = length(options_.filter_step_ahead);

MaxNumberOfBytes=options_.MaxNumberOfBytes;
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
nvobs     = size(options_.varobs,1);
iendo = 1:endo_nbr;
horizon = options_.forecast;
filtered_vars = options_.filtered_vars;
if horizon
    i_last_obs = gend+(1-M_.maximum_endo_lag:0);
end
maxlag = M_.maximum_endo_lag;

if strcmpi(type,'posterior')
    DirectoryName = CheckPath('metropolis',M_.dname);
    load_last_mh_history_file(DirectoryName, M_.fname);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine;
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); LastMhFile = TotalNumberOfMhFiles;
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    clear record;
    if ~isempty(options_.sub_draws)
        B = options_.sub_draws;
        if B > NumberOfDraws
            B = NumberOfDraws;
        end
    else
        B = min(1200, round(0.25*NumberOfDraws));
    end
elseif strcmpi(type,'gsa')
    RootDirectoryName = CheckPath('gsa',M_.dname);
    if options_.opt_gsa.pprior
        DirectoryName = CheckPath(['gsa',filesep,'prior'],M_.dname);
        load([ RootDirectoryName filesep  M_.fname '_prior.mat'],'lpmat0','lpmat','istable')
    else
        DirectoryName = CheckPath(['gsa',filesep,'mc'],M_.dname);
        load([ RootDirectoryName filesep  M_.fname '_mc.mat'],'lpmat0','lpmat','istable')
    end
    if ~isempty(lpmat0),
        x=[lpmat0(istable,:) lpmat(istable,:)];
    else
        x=lpmat(istable,:);
    end
    clear lpmat lpmat0 istable
    NumberOfDraws=size(x,1);
    B=NumberOfDraws; 
elseif strcmpi(type,'prior')
    DirectoryName = CheckPath('prior',M_.dname);
    if ~isempty(options_.subdraws)
        B = options_.subdraws;
    else
        B = 1200;
    end
end

MAX_nruns = min(B,ceil(MaxNumberOfBytes/(npar+2)/8));
MAX_nsmoo = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_ninno = min(B,ceil(MaxNumberOfBytes/(exo_nbr*gend)/8));
MAX_nerro = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8));

if naK
    MAX_naK   = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)* ...
                                             length(options_.filter_step_ahead)*gend)/8));
end

if horizon
    MAX_nforc1 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
    MAX_nforc2 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/ ...
                            8));
    IdObs    = bayestopt_.mfys;
    
end
MAX_momentsno = min(B,ceil(MaxNumberOfBytes/(get_moments_size(options_)*8)));

varlist = options_.varlist;
if isempty(varlist)
    varlist = M_.endo_names(1:M_.orig_endo_nbr, :);
end
nvar = size(varlist,1);
SelecVariables = [];
for i=1:nvar
    if ~isempty(strmatch(varlist(i,:),M_.endo_names,'exact'))
        SelecVariables = [SelecVariables;strmatch(varlist(i,:),M_.endo_names,'exact')];
    end
end

irun = ones(7,1);
ifil = zeros(7,1);


stock_param = zeros(MAX_nruns, npar);
stock_logpo = zeros(MAX_nruns,1);
stock_ys = zeros(MAX_nruns,endo_nbr);
run_smoother = 0;
if options_.smoother
    stock_smooth = zeros(endo_nbr,gend,MAX_nsmoo);
    stock_innov  = zeros(exo_nbr,gend,B);
    stock_error = zeros(nvobs,gend,MAX_nerro);
    stock_update = zeros(endo_nbr,gend,MAX_nsmoo);
    run_smoother = 1;
end

if options_.filter_step_ahead
    stock_filter_step_ahead = zeros(naK,endo_nbr,gend+ ...
                                    options_.filter_step_ahead(end),MAX_naK);
    run_smoother = 1;
end
if options_.forecast
    stock_forcst_mean = zeros(endo_nbr,horizon+maxlag,MAX_nforc1);
    stock_forcst_point = zeros(endo_nbr,horizon+maxlag,MAX_nforc2);
    run_smoother = 1;
end



% Store the variable mandatory for local/remote parallel computing.

localVars.type=type;
localVars.run_smoother=run_smoother;
localVars.gend=gend;
localVars.Y=Y;
localVars.data_index=data_index;
localVars.missing_value=missing_value;
localVars.varobs=options_.varobs;
localVars.irun=irun;
localVars.endo_nbr=endo_nbr;
localVars.nvn=nvn;
localVars.naK=naK;
localVars.horizon=horizon;
localVars.iendo=iendo;
if horizon
    localVars.i_last_obs=i_last_obs;
    localVars.IdObs=IdObs;
    localVars.MAX_nforc1=MAX_nforc1;
    localVars.MAX_nforc2=MAX_nforc2;
end
localVars.exo_nbr=exo_nbr;
localVars.maxlag=maxlag;
localVars.MAX_nsmoo=MAX_nsmoo;
localVars.MAX_ninno=MAX_ninno;
localVars.MAX_nerro = MAX_nerro;
if naK
    localVars.MAX_naK=MAX_naK;
end
localVars.MAX_nruns=MAX_nruns;
localVars.MAX_momentsno = MAX_momentsno;
localVars.ifil=ifil;
localVars.DirectoryName = DirectoryName;

if strcmpi(type,'posterior'),
    b=0;
    while b<=B
        b = b + 1;
        [x(b,:), logpost(b,1)] = GetOneDraw(type);
    end
    localVars.logpost=logpost;
end

if ~strcmpi(type,'prior'),
    localVars.x=x;
end

% Like sequential execution!
if isnumeric(options_.parallel),
    [fout] = prior_posterior_statistics_core(localVars,1,B,0);
    % Parallel execution!
else
    [nCPU, totCPU, nBlockPerCPU] = distributeJobs(options_.parallel, 1, B);
    ifil=zeros(7,totCPU);
    for j=1:totCPU-1,
        if run_smoother
            nfiles = ceil(nBlockPerCPU(j)/MAX_nsmoo);
            ifil(1,j+1) =ifil(1,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_ninno);
            ifil(2,j+1) =ifil(2,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_nerro);
            ifil(3,j+1) =ifil(3,j)+nfiles;
        end
        if naK
            nfiles = ceil(nBlockPerCPU(j)/MAX_naK);
            ifil(4,j+1) =ifil(4,j)+nfiles;
        end
        nfiles = ceil(nBlockPerCPU(j)/MAX_nruns);
        ifil(5,j+1) =ifil(5,j)+nfiles;
        if horizon
            nfiles = ceil(nBlockPerCPU(j)/MAX_nforc1);
            ifil(6,j+1) =ifil(6,j)+nfiles;
            nfiles = ceil(nBlockPerCPU(j)/MAX_nforc2);
            ifil(7,j+1) =ifil(7,j)+nfiles;
        end
    end
    localVars.ifil = ifil;
    globalVars = struct('M_',M_, ...
                        'options_', options_, ...
                        'bayestopt_', bayestopt_, ...
                        'estim_params_', estim_params_, ...
                        'oo_', oo_);    
    % which files have to be copied to run remotely
    NamFileInput(1,:) = {'',[M_.fname '_static.m']};
    NamFileInput(2,:) = {'',[M_.fname '_dynamic.m']};
    if options_.steadystate_flag,
        NamFileInput(length(NamFileInput)+1,:)={'',[M_.fname '_steadystate.m']};
    end
    [fout] = masterParallel(options_.parallel, 1, B,NamFileInput,'prior_posterior_statistics_core', localVars,globalVars, options_.parallel_info);
    
end
ifil = fout(end).ifil;

% ??????????
stock_gend=gend;
stock_data=Y;
save([DirectoryName '/' M_.fname '_data.mat'],'stock_gend','stock_data');

if strcmpi(type,'gsa')
    return
end

if ~isnumeric(options_.parallel),
    leaveSlaveOpen = options_.parallel_info.leaveSlaveOpen;
    if options_.parallel_info.leaveSlaveOpen == 0,
        % Commenting for testing!!!
        options_.parallel_info.leaveSlaveOpen = 1; % Force locally to leave open remote matlab sessions (repeated pm3 calls)
    end
end

if options_.smoother
    pm3(endo_nbr,gend,ifil(1),B,'Smoothed variables',...
        '',M_.endo_names(1:M_.orig_endo_nbr, :),M_.endo_names_tex,M_.endo_names,...
        varlist,'SmoothedVariables',DirectoryName,'_smooth');
    pm3(exo_nbr,gend,ifil(2),B,'Smoothed shocks',...
        '',M_.exo_names,M_.exo_names_tex,M_.exo_names,...
        M_.exo_names,'SmoothedShocks',DirectoryName,'_inno');
    if nvn
        % needs to  be fixed
        %        pm3(endo_nbr,gend,ifil(3),B,'Smoothed measurement errors',...
        %            M_.endo_names(SelecVariables),M_.endo_names,'tit_tex',M_.endo_names,...
        %            'names2','smooth_errors',[M_.fname '/metropolis'],'_error')
    end
end

if options_.filtered_vars
    pm3(endo_nbr,gend,ifil(1),B,'Updated Variables',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'UpdatedVariables',DirectoryName, ...
        '_update');
    pm3(endo_nbr,gend,ifil(4),B,'One step ahead forecast (filtered variables)',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'FilteredVariables',DirectoryName,'_filter_step_ahead');
end

if options_.forecast
    pm3(endo_nbr,horizon+maxlag,ifil(6),B,'Forecasted variables (mean)',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'MeanForecast',DirectoryName,'_forc_mean');
    pm3(endo_nbr,horizon+maxlag,ifil(6),B,'Forecasted variables (point)',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'PointForecast',DirectoryName,'_forc_point');
end


if ~isnumeric(options_.parallel),
    options_.parallel_info.leaveSlaveOpen = leaveSlaveOpen;
    if leaveSlaveOpen == 0,
        closeSlave(options_.parallel,options_.parallel_info.RemoteTmpFolder),
    end
end
