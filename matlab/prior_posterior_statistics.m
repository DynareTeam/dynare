function prior_posterior_statistics(type,dataset,dataset_info)

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
% 
% PARALLEL CONTEXT
% See the comments in the posterior_sampler.m funtion.
% 
% 
% Copyright (C) 2005-2016 Dynare Team
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

Y = transpose(dataset.data);
gend = dataset.nobs;
data_index = dataset_info.missing.aindex;
missing_value = dataset_info.missing.state;
mean_varobs = dataset_info.descriptive.mean;


nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
naK = length(options_.filter_step_ahead);

MaxNumberOfBytes=options_.MaxNumberOfBytes;
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
meas_err_nbr=length(M_.Correlation_matrix_ME);
iendo = 1:endo_nbr;
horizon = options_.forecast;
IdObs    = bayestopt_.mfys;
if horizon
    i_last_obs = gend+(1-M_.maximum_endo_lag:0);
end
maxlag = M_.maximum_endo_lag;

if strcmpi(type,'posterior')
    DirectoryName = CheckPath('metropolis',M_.dname);
    B = options_.sub_draws;
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
    B = size(x,1);
elseif strcmpi(type,'prior')
    DirectoryName = CheckPath('prior',M_.dname);
    B = options_.prior_draws;
end

MAX_nruns = min(B,ceil(MaxNumberOfBytes/(npar+2)/8));
MAX_nsmoo = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_n_smoothed_constant = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_n_smoothed_trend = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_n_trend_coeff = min(B,ceil(MaxNumberOfBytes/endo_nbr/8));
MAX_ninno = min(B,ceil(MaxNumberOfBytes/(exo_nbr*gend)/8));
MAX_nerro = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8));

if naK
    MAX_naK   = min(B,ceil(MaxNumberOfBytes/(endo_nbr* ...
                                             length(options_.filter_step_ahead)*(gend+max(options_.filter_step_ahead)))/8));
end

if horizon
    MAX_nforc1 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
    MAX_nforc2 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
    if ~isequal(M_.H,0)
        MAX_nforc_ME = min(B,ceil(MaxNumberOfBytes/((size(options_.varobs,1))*(horizon+maxlag))/8));
    end
end
MAX_momentsno = min(B,ceil(MaxNumberOfBytes/(get_moments_size(options_)*8)));

if options_.filter_covariance
    MAX_filter_covariance = min(B,ceil(MaxNumberOfBytes/(endo_nbr^2*(gend+1))/8));
end

if options_.smoothed_state_uncertainty 
    MAX_n_smoothed_state_uncertainty = min(B,ceil(MaxNumberOfBytes/((endo_nbr*endo_nbr)*gend)/8));
end

varlist = options_.varlist;
if isempty(varlist)
    varlist = char(sort(cellstr(M_.endo_names(1:M_.orig_endo_nbr,:))));
end
nvar = size(varlist,1);
SelecVariables = [];
for i=1:nvar
    if ~isempty(strmatch(varlist(i,:),M_.endo_names,'exact'))
        SelecVariables = [SelecVariables;strmatch(varlist(i,:),M_.endo_names,'exact')];
    end
end

n_variables_to_fill=13;

irun = ones(n_variables_to_fill,1);
ifil = zeros(n_variables_to_fill,1);

run_smoother = 0;
if options_.smoother || options_.forecast || options_.filter_step_ahead || options_.smoothed_state_uncertainty
    run_smoother = 1;
end

filter_covariance=0;
if options_.filter_covariance
    filter_covariance=1;
end

smoothed_state_uncertainty=0;
if options_.smoothed_state_uncertainty
    smoothed_state_uncertainty=1;
end

% Store the variable mandatory for local/remote parallel computing.

localVars.type=type;
localVars.run_smoother=run_smoother;
localVars.filter_covariance=filter_covariance;
localVars.smoothed_state_uncertainty=smoothed_state_uncertainty;
localVars.gend=gend;
localVars.Y=Y;
localVars.data_index=data_index;
localVars.missing_value=missing_value;
localVars.varobs=options_.varobs;
localVars.mean_varobs=mean_varobs;
localVars.irun=irun;
localVars.endo_nbr=endo_nbr;
localVars.nvn=nvn;
localVars.naK=naK;
localVars.horizon=horizon;
localVars.iendo=iendo;
localVars.IdObs=IdObs;
if horizon
    localVars.i_last_obs=i_last_obs;
    localVars.MAX_nforc1=MAX_nforc1;
    localVars.MAX_nforc2=MAX_nforc2;
    if ~isequal(M_.H,0)
        localVars.MAX_nforc_ME = MAX_nforc_ME;
    end
end
localVars.exo_nbr=exo_nbr;
localVars.maxlag=maxlag;
localVars.MAX_nsmoo=MAX_nsmoo;
localVars.MAX_ninno=MAX_ninno;
localVars.MAX_nerro = MAX_nerro;
if naK
    localVars.MAX_naK=MAX_naK;
end
if options_.filter_covariance
    localVars.MAX_filter_covariance = MAX_filter_covariance;
end
if options_.smoothed_state_uncertainty 
    localVars.MAX_n_smoothed_state_uncertainty = MAX_n_smoothed_state_uncertainty ;
end
localVars.MAX_n_smoothed_constant=MAX_n_smoothed_constant;
localVars.MAX_n_smoothed_trend=MAX_n_smoothed_trend;
localVars.MAX_n_trend_coeff=MAX_n_trend_coeff;
localVars.MAX_nruns=MAX_nruns;
localVars.MAX_momentsno = MAX_momentsno;
localVars.ifil=ifil;
localVars.DirectoryName = DirectoryName;

if strcmpi(type,'posterior'),
    BaseName = [DirectoryName filesep M_.fname];
    load_last_mh_history_file(DirectoryName, M_.fname);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine; 
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); 
    LastMhFile = TotalNumberOfMhFiles; 
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    mh_nblck = options_.mh_nblck;
    if B==NumberOfDraws*mh_nblck,
        % we load all retained MH runs !
        logpost=GetAllPosteriorDraws(0, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws);
        for column=1:npar
            x(:,column) = GetAllPosteriorDraws(column, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws);
        end
    else

    b=0;
    logpost=NaN(B,1);
        while b<B
        b = b + 1;
            [xtmp, ltmp] = GetOneDraw(type);
            if ~ismember(ltmp,logpost(1:b-1))
                x(b,:)=xtmp;
                logpost(b,1) = ltmp;
            else
                b=b-1;
            end
        end
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
    ifil=zeros(n_variables_to_fill,totCPU);
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
            if ~isequal(M_.H,0)
                nfiles = ceil(nBlockPerCPU(j)/MAX_nforc_ME);
                ifil(12,j+1) =ifil(12,j)+nfiles;
            end
        end
        if options_.filter_covariance
            nfiles = ceil(nBlockPerCPU(j)/MAX_filter_covariance);
            ifil(8,j+1) =ifil(8,j)+nfiles;
        end
        if run_smoother
            nfiles = ceil(nBlockPerCPU(j)/MAX_n_trend_coeff);
            ifil(9,j+1) =ifil(9,j)+nfiles;  
            nfiles = ceil(nBlockPerCPU(j)/MAX_n_smoothed_constant);
            ifil(10,j+1) =ifil(10,j)+nfiles;  
            nfiles = ceil(nBlockPerCPU(j)/MAX_n_smoothed_trend);
            ifil(11,j+1) =ifil(11,j)+nfiles;  
            if smoothed_state_uncertainty
                nfiles = ceil(nBlockPerCPU(j)/MAX_n_smoothed_state_uncertainty);
                ifil(13,j+1) =ifil(13,j)+nfiles;           
            end
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
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'SmoothedVariables',DirectoryName,'_smooth');
    pm3(exo_nbr,gend,ifil(2),B,'Smoothed shocks',...
        '',M_.exo_names,M_.exo_names_tex,M_.exo_names,...
        M_.exo_names,'SmoothedShocks',DirectoryName,'_inno');
    pm3(endo_nbr,1,ifil(9),B,'Trend_coefficients',...
    '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'TrendCoeff',DirectoryName,'_trend_coeff');
    pm3(endo_nbr,gend,ifil(10),B,'Smoothed constant',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'Constant',DirectoryName,'_smoothed_constant');
    pm3(endo_nbr,gend,ifil(11),B,'Smoothed trend',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'Trend',DirectoryName,'_smoothed_trend');
    pm3(endo_nbr,gend,ifil(1),B,'Updated Variables',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'UpdatedVariables',DirectoryName, ...
        '_update');
    if smoothed_state_uncertainty
        pm3(endo_nbr,endo_nbr,ifil(13),B,'State Uncertainty',...
            '',varlist,M_.endo_names_tex,M_.endo_names,...
            varlist,'StateUncertainty',DirectoryName,'_state_uncert');
    end
    
    if nvn
        for obs_iter=1:length(options_.varobs)        
            meas_error_names{obs_iter,1}=['SE_EOBS_' M_.endo_names(strmatch(options_.varobs{obs_iter},M_.endo_names,'exact'),:)];
            texnames{obs_iter,1}=['SE_EOBS_' M_.endo_names_tex(strmatch(options_.varobs{obs_iter},M_.endo_names,'exact'),:)];
        end
        meas_error_names=char(meas_error_names);
        texnames=char(texnames);
        pm3(meas_err_nbr,gend,ifil(3),B,'Smoothed measurement errors',...
           '',meas_error_names,texnames,meas_error_names,...
           meas_error_names,'SmoothedMeasurementErrors',DirectoryName,'_error')
    end
end

if options_.filtered_vars
    pm3(endo_nbr,gend,ifil(4),B,'One step ahead forecast (filtered variables)',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'FilteredVariables',DirectoryName,'_filter_step_ahead');
end

if options_.forecast
    pm3(endo_nbr,horizon,ifil(6),B,'Forecasted variables (mean)',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'MeanForecast',DirectoryName,'_forc_mean');
    pm3(endo_nbr,horizon,ifil(7),B,'Forecasted variables (point)',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'PointForecast',DirectoryName,'_forc_point');
    if ~isequal(M_.H,0) && ~isempty(intersect(options_.varobs,varlist))
        texnames=[];
        for obs_iter=1:length(options_.varobs)        
            obs_names{obs_iter,1}=M_.endo_names(strmatch(options_.varobs{obs_iter},M_.endo_names,'exact'),:);
            texnames{obs_iter,1}=M_.endo_names_tex(strmatch(options_.varobs{obs_iter},M_.endo_names,'exact'),:);
        end
        obs_names=char(obs_names);
        texnames=char(texnames);
        varlist_forecast_ME=intersect(options_.varobs,varlist);
        pm3(meas_err_nbr,horizon,ifil(12),B,'Forecasted variables (point) with ME',...
           '',char(varlist_forecast_ME),texnames,obs_names,...
           char(varlist_forecast_ME),'PointForecastME',DirectoryName,'_forc_point_ME')    
    end
end

if options_.filter_covariance
    pm3(endo_nbr,endo_nbr,ifil(8),B,'Filtered covariances',...
        '',varlist,M_.endo_names_tex,M_.endo_names,...
        varlist,'FilterCovariance',DirectoryName,'_filter_covar');
end


if ~isnumeric(options_.parallel),
    options_.parallel_info.leaveSlaveOpen = leaveSlaveOpen;
    if leaveSlaveOpen == 0,
        closeSlave(options_.parallel,options_.parallel_info.RemoteTmpFolder),
    end
end