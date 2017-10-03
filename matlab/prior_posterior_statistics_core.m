function myoutput=prior_posterior_statistics_core(myinputs,fpar,B,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% Core functionality for prior_posterior.m function, which can be parallelized.
% See also the comment in posterior_sampler_core.m funtion.
%
% INPUTS
%   See the comment in posterior_sampler_core.m funtion.
%
% OUTPUTS
% o myoutput  [struc]
%  Contained OutputFileName_smooth;
%                          _update;
%                          _inno;
%                          _error;
%                          _filter_step_ahead;
%                          _param;
%                          _forc_mean;
%                          _forc_point;
%                          _forc_point_ME;
%                          _filter_covar;
%                          _trend_coeff;
%                          _smoothed_trend;
%                          _smoothed_constant;
%                          _state_uncert;
%
% ALGORITHM
%   Portion of prior_posterior.m function.
% This file is part of Dynare.
%
% SPECIAL REQUIREMENTS.
%   None.

% Copyright (C) 2005-2017 Dynare Team
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

global options_ oo_ M_ bayestopt_ estim_params_

if nargin<4
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

type=myinputs.type;
run_smoother=myinputs.run_smoother;
filter_covariance=myinputs.filter_covariance;
smoothed_state_uncertainty=myinputs.smoothed_state_uncertainty;
gend=myinputs.gend;
Y=myinputs.Y;
data_index=myinputs.data_index;
missing_value=myinputs.missing_value;
varobs=myinputs.varobs;
mean_varobs=myinputs.mean_varobs;
irun=myinputs.irun;
endo_nbr=myinputs.endo_nbr;
nvn=myinputs.nvn;
naK=myinputs.naK;
horizon=myinputs.horizon;
iendo=myinputs.iendo;
IdObs=myinputs.IdObs; %index of observables
if horizon
    i_last_obs=myinputs.i_last_obs;
    MAX_nforc1=myinputs.MAX_nforc1;
    MAX_nforc2=myinputs.MAX_nforc2;
    if ~isequal(M_.H,0)
        MAX_nforc_ME=myinputs.MAX_nforc_ME;
    end
end
if naK
    MAX_naK=myinputs.MAX_naK;
end
if filter_covariance
    MAX_filter_covariance=myinputs.MAX_filter_covariance;
end
if smoothed_state_uncertainty
    MAX_n_smoothed_state_uncertainty=myinputs.MAX_n_smoothed_state_uncertainty;
end

exo_nbr=myinputs.exo_nbr;
maxlag=myinputs.maxlag;
MAX_nsmoo=myinputs.MAX_nsmoo;
MAX_ninno=myinputs.MAX_ninno;
MAX_n_smoothed_constant=myinputs.MAX_n_smoothed_constant;
MAX_n_smoothed_trend=myinputs.MAX_n_smoothed_trend;
MAX_n_trend_coeff=myinputs.MAX_n_trend_coeff;
MAX_nerro = myinputs.MAX_nerro;
MAX_nruns=myinputs.MAX_nruns;
MAX_momentsno = myinputs.MAX_momentsno;
ifil=myinputs.ifil;

if ~strcmpi(type,'prior')
    x=myinputs.x;
    if strcmpi(type,'posterior')
        logpost=myinputs.logpost;
    end
end
if whoiam
    Parallel=myinputs.Parallel;
end

% DirectoryName = myinputs.DirectoryName;
if strcmpi(type,'posterior')
    DirectoryName = CheckPath('metropolis',M_.dname);
elseif strcmpi(type,'gsa')
    if options_.opt_gsa.pprior
        DirectoryName = CheckPath(['gsa',filesep,'prior'],M_.dname);
    else
        DirectoryName = CheckPath(['gsa',filesep,'mc'],M_.dname);
    end
elseif strcmpi(type,'prior')
    DirectoryName = CheckPath('prior',M_.dname);
end

RemoteFlag = 0;
if whoiam
    if Parallel(ThisMatlab).Local==0
        RemoteFlag =1;
    end
    ifil=ifil(:,whoiam);
    prct0={0,whoiam,Parallel(ThisMatlab)};
else
    prct0=0;
end
h = dyn_waitbar(prct0,['Taking ',type,' subdraws...']);

if RemoteFlag==1
    OutputFileName_smooth = {};
    OutputFileName_update = {};
    OutputFileName_inno = {};
    OutputFileName_error = {};
    OutputFileName_filter_step_ahead = {};
    OutputFileName_param = {};
    OutputFileName_forc_mean = {};
    OutputFileName_forc_point = {};
    OutputFileName_forc_point_ME = {};
    OutputFileName_filter_covar ={};
    OutputFileName_trend_coeff = {};
    OutputFileName_smoothed_trend = {};
    OutputFileName_smoothed_constant = {};
    % OutputFileName_moments = {};
end

%initialize arrays
if run_smoother
    stock_smooth=NaN(endo_nbr,gend,MAX_nsmoo);
    stock_update=NaN(endo_nbr,gend,MAX_nsmoo);
    stock_innov=NaN(M_.exo_nbr,gend,MAX_ninno);
    stock_smoothed_constant=NaN(endo_nbr,gend,MAX_n_smoothed_constant);
    stock_smoothed_trend=NaN(endo_nbr,gend,MAX_n_smoothed_trend);
    stock_trend_coeff = zeros(endo_nbr,MAX_n_trend_coeff);
    if horizon
        stock_forcst_mean= NaN(endo_nbr,horizon,MAX_nforc1);
        stock_forcst_point = NaN(endo_nbr,horizon,MAX_nforc2);
        if ~isequal(M_.H,0)
            stock_forcst_point_ME = NaN(length(varobs),horizon,MAX_nforc_ME);
        end
    end
end
if nvn
    stock_error = NaN(length(varobs),gend,MAX_nerro);
end
if naK
    stock_filter_step_ahead =NaN(length(options_.filter_step_ahead),endo_nbr,gend+max(options_.filter_step_ahead),MAX_naK);
end
stock_param = NaN(MAX_nruns,size(myinputs.x,2));
stock_logpo = NaN(MAX_nruns,1);
stock_ys = NaN(MAX_nruns,endo_nbr);
if filter_covariance
    stock_filter_covariance = zeros(endo_nbr,endo_nbr,gend+1,MAX_filter_covariance);
end
if smoothed_state_uncertainty
    stock_smoothed_uncert = zeros(endo_nbr,endo_nbr,gend,MAX_n_smoothed_state_uncertainty);
end

for b=fpar:B
    if strcmpi(type,'prior')

        [deep, logpo] = GetOneDraw(type,M_,estim_params_,oo_,options_,bayestopt_);

    else
        deep = x(b,:);
        if strcmpi(type,'posterior')
            logpo = logpost(b);
        else
            logpo = evaluate_posterior_kernel(deep',M_,estim_params_,oo_,options_,bayestopt_);
        end
    end
    M_ = set_all_parameters(deep,estim_params_,M_);

    if run_smoother
        [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        [alphahat,etahat,epsilonhat,alphatilde,SteadyState,trend_coeff,aK,junk1,junk2,P,junk4,junk5,trend_addition,state_uncertainty,M_,oo_,options_,bayestopt_] = ...
            DsgeSmoother(deep,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);

        stock_trend_coeff(options_.varobs_id,irun(9))=trend_coeff;
        stock_smoothed_trend(IdObs,:,irun(11))=trend_addition;
        if options_.loglinear %reads values from smoother results, which are in dr-order and put them into declaration order
            constant_part=repmat(log(SteadyState(dr.order_var)),1,gend);
            stock_smooth(dr.order_var,:,irun(1)) = alphahat(1:endo_nbr,:)+ ...
                constant_part;
            stock_update(dr.order_var,:,irun(1)) = alphatilde(1:endo_nbr,:)+ ...
                constant_part;
        else
            constant_part=repmat(SteadyState(dr.order_var),1,gend);
            stock_smooth(dr.order_var,:,irun(1)) = alphahat(1:endo_nbr,:)+ ...
                constant_part;
            stock_update(dr.order_var,:,irun(1)) = alphatilde(1:endo_nbr,:)+ ...
                constant_part;
        end
        stock_smoothed_constant(dr.order_var,:,irun(10))=constant_part;
        %% Compute constant for observables
        if options_.prefilter == 1 %as mean is taken after log transformation, no distinction is needed here
            constant_part=repmat(mean_varobs',1,gend);
        elseif options_.prefilter == 0 && options_.loglinear == 1 %logged steady state must be used
            constant_part=repmat(log(SteadyState(IdObs)),1,gend);
        elseif options_.prefilter == 0 && options_.loglinear == 0 %unlogged steady state must be used
            constant_part=repmat(SteadyState(IdObs),1,gend);
        end
        %add trend to observables
        if options_.prefilter
            %do correction for prefiltering for observed variables
            if options_.loglinear
                mean_correction=-repmat(log(SteadyState(IdObs)),1,gend)+constant_part;
            else
                mean_correction=-repmat(SteadyState(IdObs),1,gend)+constant_part;
            end
            stock_smoothed_constant(IdObs,:,irun(10))=stock_smoothed_constant(IdObs,:,irun(10))+mean_correction;
            %smoothed variables are E_T(y_t) so no trend shift is required
            stock_smooth(IdObs,:,irun(1))=stock_smooth(IdObs,:,irun(1))+trend_addition+mean_correction;
            %updated variables are E_t(y_t) so no trend shift is required
            stock_update(IdObs,:,irun(1))=stock_update(IdObs,:,irun(1))+trend_addition+mean_correction;
        else
            stock_smooth(IdObs,:,irun(1))=stock_smooth(IdObs,:,irun(1))+trend_addition;
            stock_update(IdObs,:,irun(1))=stock_update(IdObs,:,irun(1))+trend_addition;
        end
        stock_innov(:,:,irun(2))  = etahat;
        if nvn
            stock_error(:,:,irun(3))  = epsilonhat;
        end
        if naK
            %filtered variable E_t(y_t+k) requires to shift trend by k periods
            %write variables into declaration order
            if options_.loglinear %reads values from smoother results, which are in dr-order and put them into declaration order
                constant_part=repmat(log(SteadyState(dr.order_var))',[length(options_.filter_step_ahead),1,gend+max(options_.filter_step_ahead)]);
            else
                constant_part=repmat(SteadyState(dr.order_var)',[length(options_.filter_step_ahead),1,gend+max(options_.filter_step_ahead)]);
            end
            stock_filter_step_ahead(:,dr.order_var,:,irun(4)) = aK(options_.filter_step_ahead,1:endo_nbr,:) + constant_part;
            %now add trend to observables
            for ii=1:length(options_.filter_step_ahead)
                if options_.prefilter
                    zdim = size(stock_filter_step_ahead(ii,IdObs,:,irun(4)));
                    squeezed = reshape(stock_filter_step_ahead(ii,IdObs,:,irun(4)), [zdim(2:end) 1]);
                    stock_filter_step_ahead(ii,IdObs,:,irun(4)) = squeezed ...
                        +repmat(mean_correction(:,1),1,gend+max(options_.filter_step_ahead)) ... %constant correction
                        +[trend_addition repmat(trend_addition(:,end),1,max(options_.filter_step_ahead))+trend_coeff*[1:max(options_.filter_step_ahead)]]; %trend
                else
                    zdim = size(stock_filter_step_ahead(ii,IdObs,:,irun(4)));
                    squeezed = reshape(stock_filter_step_ahead(ii,IdObs,:,irun(4)), [zdim(2:end) 1]);
                    stock_filter_step_ahead(ii,IdObs,:,irun(4)) = squeezed ...
                        +[trend_addition repmat(trend_addition(:,end),1,max(options_.filter_step_ahead))+trend_coeff*[1:max(options_.filter_step_ahead)]]; %trend
                end
            end
        end
        if horizon
            yyyy = alphahat(iendo,i_last_obs);
            yf = forcst2a(yyyy,dr,zeros(horizon,exo_nbr));
            if options_.prefilter
                % add mean
                yf(:,IdObs) = yf(:,IdObs)+repmat(mean_varobs, ...
                                                 horizon+maxlag,1);
                % add trend, taking into account that last point of sample is still included in forecasts and only cut off later
                yf(:,IdObs) = yf(:,IdObs)+((options_.first_obs-1)+gend+[1-maxlag:horizon]')*trend_coeff'-...
                    repmat(mean(trend_coeff*[options_.first_obs:options_.first_obs+gend-1],2)',length(1-maxlag:horizon),1); %center trend
            else
                % add trend, taking into account that last point of sample is still included in forecasts and only cut off later
                yf(:,IdObs) = yf(:,IdObs)+((options_.first_obs-1)+gend+[1-maxlag:horizon]')*trend_coeff';
            end
            if options_.loglinear
                yf = yf+repmat(log(SteadyState'),horizon+maxlag,1);
            else
                yf = yf+repmat(SteadyState',horizon+maxlag,1);
            end
            yf1 = forcst2(yyyy,horizon,dr,1);
            if options_.prefilter == 1
                % add mean
                yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
                    repmat(mean_varobs,[horizon+maxlag,1,1]);
                % add trend, taking into account that last point of sample is still included in forecasts and only cut off later
                yf1(:,IdObs) = yf1(:,IdObs)+((options_.first_obs-1)+gend+[1-maxlag:horizon]')*trend_coeff'-...
                    repmat(mean(trend_coeff*[options_.first_obs:options_.first_obs+gend-1],2)',length(1-maxlag:horizon),1); %center trend
            else
                % add trend, taking into account that last point of sample is still included in forecasts and only cut off later
                yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat(((options_.first_obs-1)+gend+[1-maxlag:horizon]')* ...
                                                       trend_coeff',[1,1,1]);
            end
            if options_.loglinear
                yf1 = yf1 + repmat(log(SteadyState'),[horizon+maxlag,1,1]);
            else
                yf1 = yf1 + repmat(SteadyState',[horizon+maxlag,1,1]);
            end

            stock_forcst_mean(:,:,irun(6)) = yf(maxlag+1:end,:)';
            stock_forcst_point(:,:,irun(7)) = yf1(maxlag+1:end,:)';
            if ~isequal(M_.H,0)
                ME_shocks=zeros(length(varobs),horizon);
                i_exo_var = setdiff([1:length(varobs)],find(diag(M_.H) == 0));
                nxs = length(i_exo_var);
                chol_H = chol(M_.H(i_exo_var,i_exo_var));
                if ~isempty(M_.H)
                    ME_shocks(i_exo_var,:) = chol_H*randn(nxs,horizon);
                end
                stock_forcst_point_ME(:,:,irun(12)) = yf1(maxlag+1:end,IdObs)'+ME_shocks;
            end
        end
        if filter_covariance
            stock_filter_covariance(dr.order_var,dr.order_var,:,irun(8)) = P;
        end
        if smoothed_state_uncertainty
            stock_smoothed_uncert(dr.order_var,dr.order_var,:,irun(13)) = state_uncertainty;
        end
    else
        [T,R,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);
    end
    stock_param(irun(5),:) = deep;
    stock_logpo(irun(5),1) = logpo;
    stock_ys(irun(5),:) = SteadyState';


    irun = irun +  ones(13,1);


    if run_smoother && (irun(1) > MAX_nsmoo || b == B)
        stock = stock_smooth(:,:,1:irun(1)-1);
        ifil(1) = ifil(1) + 1;
        save([DirectoryName '/' M_.fname '_smooth' int2str(ifil(1)) '.mat'],'stock');

        stock = stock_update(:,:,1:irun(1)-1);
        save([DirectoryName '/' M_.fname '_update' int2str(ifil(1)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_smooth = [OutputFileName_smooth; {[DirectoryName filesep], [M_.fname '_smooth' int2str(ifil(1)) '.mat']}];
            OutputFileName_update = [OutputFileName_update; {[DirectoryName filesep], [M_.fname '_update' int2str(ifil(1)) '.mat']}];
        end
        irun(1) = 1;
    end

    if run_smoother && (irun(2) > MAX_ninno || b == B)
        stock = stock_innov(:,:,1:irun(2)-1);
        ifil(2) = ifil(2) + 1;
        save([DirectoryName '/' M_.fname '_inno' int2str(ifil(2)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_inno = [OutputFileName_inno; {[DirectoryName filesep], [M_.fname '_inno' int2str(ifil(2)) '.mat']}];
        end
        irun(2) = 1;
    end

    if run_smoother && nvn && (irun(3) > MAX_nerro || b == B)
        stock = stock_error(:,:,1:irun(3)-1);
        ifil(3) = ifil(3) + 1;
        save([DirectoryName '/' M_.fname '_error' int2str(ifil(3)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_error = [OutputFileName_error; {[DirectoryName filesep], [M_.fname '_error' int2str(ifil(3)) '.mat']}];
        end
        irun(3) = 1;
    end

    if run_smoother && naK && (irun(4) > MAX_naK || b == B)
        stock = stock_filter_step_ahead(:,:,:,1:irun(4)-1);
        ifil(4) = ifil(4) + 1;
        save([DirectoryName '/' M_.fname '_filter_step_ahead' int2str(ifil(4)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_filter_step_ahead = [OutputFileName_filter_step_ahead; {[DirectoryName filesep], [M_.fname '_filter_step_ahead' int2str(ifil(4)) '.mat']}];
        end
        irun(4) = 1;
    end

    if irun(5) > MAX_nruns || b == B
        stock = stock_param(1:irun(5)-1,:);
        stock_logpo = stock_logpo(1:irun(5)-1);
        stock_ys = stock_ys(1:irun(5)-1,:);
        ifil(5) = ifil(5) + 1;
        save([DirectoryName '/' M_.fname '_param' int2str(ifil(5)) '.mat'],'stock','stock_logpo','stock_ys');
        if RemoteFlag==1
            OutputFileName_param = [OutputFileName_param; {[DirectoryName filesep], [M_.fname '_param' int2str(ifil(5)) '.mat']}];
        end
        irun(5) = 1;
    end

    if run_smoother && horizon && (irun(6) > MAX_nforc1 || b == B)
        stock = stock_forcst_mean(:,:,1:irun(6)-1);
        ifil(6) = ifil(6) + 1;
        save([DirectoryName '/' M_.fname '_forc_mean' int2str(ifil(6)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_forc_mean = [OutputFileName_forc_mean; {[DirectoryName filesep], [M_.fname '_forc_mean' int2str(ifil(6)) '.mat']}];
        end
        irun(6) = 1;
    end

    if run_smoother && horizon && (irun(7) > MAX_nforc2 ||  b == B)
        stock = stock_forcst_point(:,:,1:irun(7)-1);
        ifil(7) = ifil(7) + 1;
        save([DirectoryName '/' M_.fname '_forc_point' int2str(ifil(7)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_forc_point = [OutputFileName_forc_point; {[DirectoryName filesep], [M_.fname '_forc_point' int2str(ifil(7)) '.mat']}];
        end
        irun(7) = 1;
    end

    if run_smoother && filter_covariance && (irun(8) > MAX_filter_covariance || b == B)
        stock = stock_filter_covariance(:,:,:,1:irun(8)-1);
        ifil(8) = ifil(8) + 1;
        save([DirectoryName '/' M_.fname '_filter_covar' int2str(ifil(8)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_filter_covar = [OutputFileName_filter_covar; {[DirectoryName filesep], [M_.fname '_filter_covar' int2str(ifil(8)) '.mat']}];
        end
        irun(8) = 1;
    end

    irun_index=9;
    if run_smoother && (irun(irun_index) > MAX_n_trend_coeff || b == B)
        stock = stock_trend_coeff(:,1:irun(irun_index)-1);
        ifil(irun_index) = ifil(irun_index) + 1;
        save([DirectoryName '/' M_.fname '_trend_coeff' int2str(ifil(irun_index)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_trend_coeff = [OutputFileName_trend_coeff; {[DirectoryName filesep], [M_.fname '_trend_coeff' int2str(ifil(irun_index)) '.mat']}];
        end
        irun(irun_index) = 1;
    end

    irun_index=10;
    if run_smoother && (irun(irun_index) > MAX_n_smoothed_constant || b == B)
        stock = stock_smoothed_constant(:,:,1:irun(irun_index)-1);
        ifil(irun_index) = ifil(irun_index) + 1;
        save([DirectoryName '/' M_.fname '_smoothed_constant' int2str(ifil(irun_index)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_smoothed_constant = [OutputFileName_smoothed_constant; {[DirectoryName filesep], [M_.fname '_smoothed_constant' int2str(ifil(irun_index)) '.mat']}];
        end
        irun(irun_index) = 1;
    end

    irun_index=11;
    if run_smoother && (irun(irun_index) > MAX_n_smoothed_trend || b == B)
        stock = stock_smoothed_trend(:,:,1:irun(irun_index)-1);
        ifil(irun_index) = ifil(irun_index) + 1;
        save([DirectoryName '/' M_.fname '_smoothed_trend' int2str(ifil(irun_index)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_smoothed_trend = [OutputFileName_smoothed_trend; {[DirectoryName filesep], [M_.fname '_smoothed_trend' int2str(ifil(irun_index)) '.mat']}];
        end
        irun(irun_index) = 1;
    end

    irun_index=12;
    if run_smoother && horizon && ~isequal(M_.H,0) && (irun(irun_index) > MAX_nforc_ME ||  b == B)
        stock = stock_forcst_point_ME(:,:,1:irun(irun_index)-1);
        ifil(irun_index) = ifil(irun_index) + 1;
        save([DirectoryName '/' M_.fname '_forc_point_ME' int2str(ifil(irun_index)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_forc_point_ME = [OutputFileName_forc_point_ME; {[DirectoryName filesep], [M_.fname '_forc_point_ME' int2str(ifil(irun_index)) '.mat']}];
        end
        irun(irun_index) = 1;
    end

    irun_index=13;
    if run_smoother && smoothed_state_uncertainty && (irun(irun_index) > MAX_n_smoothed_state_uncertainty || b == B)
        stock = stock_smoothed_uncert(:,:,:,1:irun(irun_index)-1);
        ifil(irun_index) = ifil(irun_index) + 1;
        save([DirectoryName '/' M_.fname '_state_uncert' int2str(ifil(irun_index)) '.mat'],'stock');
        if RemoteFlag==1
            OutputFileName_state_uncert = [OutputFileName_state_uncert; {[DirectoryName filesep], [M_.fname '_state_uncert' int2str(ifil(irun_index)) '.mat']}];
        end
        irun(irun_index) = 1;
    end

    dyn_waitbar((b-fpar+1)/(B-fpar+1),h);
end

myoutput.ifil=ifil;
if RemoteFlag==1
    myoutput.OutputFileName = [OutputFileName_smooth;
                        OutputFileName_update;
                        OutputFileName_inno;
                        OutputFileName_error;
                        OutputFileName_filter_step_ahead;
                        OutputFileName_param;
                        OutputFileName_forc_mean;
                        OutputFileName_forc_point;
                        OutputFileName_forc_point_ME;
                        OutputFileName_filter_covar;
                        OutputFileName_trend_coeff;
                        OutputFileName_smoothed_trend;
                        OutputFileName_smoothed_constant;
                        OutputFileName_state_uncert];
end

dyn_waitbar_close(h);
