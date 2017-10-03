function imcforecast(constrained_paths, constrained_vars, options_cond_fcst)
% Computes conditional forecasts.
%
% INPUTS
%  o constrained_paths    [double]      m*p array, where m is the number of constrained endogenous variables and p is the number of constrained periods.
%  o constrained_vars     [char]        m*x array holding the names of the controlled endogenous variables.
%  o options_cond_fcst    [structure]   containing the options. The fields are:
%                                                             + replic              [integer]   scalar, number of monte carlo simulations.
%                                                             + parameter_set       [char]      values of the estimated parameters:
%                                                                                               "posterior_mode",
%                                                                                               "posterior_mean",
%                                                                                               "posterior_median",
%                                                                                               "prior_mode" or
%                                                                                               "prior mean".
%                                                                                   [double]     np*1 array, values of the estimated parameters.
%                                                             + controlled_varexo   [char]       m*x array, list of controlled exogenous variables.
%                                                             + conf_sig            [double]     scalar in [0,1], probability mass covered by the confidence bands.
%
% OUTPUTS
%  None.
%
% SPECIAL REQUIREMENTS
%  This routine has to be called after an estimation statement or an estimated_params block.
%
% REMARKS
%  [1] Results are stored in a structure which is saved in a mat file called conditional_forecasts.mat.
%  [2] Use the function plot_icforecast to plot the results.

% Copyright (C) 2006-2017 Dynare Team
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

if ~isfield(options_cond_fcst,'parameter_set') || isempty(options_cond_fcst.parameter_set)
    if isfield(oo_,'posterior_mode')
        options_cond_fcst.parameter_set = 'posterior_mode';
    elseif isfield(oo_,'mle_mode')
        options_cond_fcst.parameter_set = 'mle_mode';
    else
        error('No valid parameter set found')
    end
end

if ~isfield(options_cond_fcst,'replic') || isempty(options_cond_fcst.replic)
    options_cond_fcst.replic = 5000;
end

if ~isfield(options_cond_fcst,'periods') || isempty(options_cond_fcst.periods)
    options_cond_fcst.periods = 40;
end

if ~isfield(options_cond_fcst,'conditional_forecast') || ~isfield(options_cond_fcst.conditional_forecast,'conf_sig')  || isempty(options_cond_fcst.conditional_forecast.conf_sig)
    options_cond_fcst.conditional_forecast.conf_sig = .8;
end

if isequal(options_cond_fcst.parameter_set,'calibration')
    estimated_model = 0;
else
    estimated_model = 1;
end

if estimated_model
    if options_.prefilter
        error('imcforecast:: Conditional forecasting does not support the prefiltering option')
    end
    if ischar(options_cond_fcst.parameter_set)
        switch options_cond_fcst.parameter_set
          case 'posterior_mode'
            xparam = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
            graph_title='Posterior Mode';
          case 'posterior_mean'
            xparam = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
            graph_title='Posterior Mean';
          case 'posterior_median'
            xparam = get_posterior_parameters('median',M_,estim_params_,oo_,options_);
            graph_title='Posterior Median';
          case 'mle_mode'
            xparam = get_posterior_parameters('mode',M_,estim_params_,oo_,options_,'mle_');
            graph_title='ML Mode';
          case 'prior_mode'
            xparam = bayestopt_.p5(:);
            graph_title='Prior Mode';
          case 'prior_mean'
            xparam = bayestopt_.p1;
            graph_title='Prior Mean';
          otherwise
            disp('imcforecast:: If the input argument is a string, then it has to be equal to:')
            disp('                   ''calibration'', ')
            disp('                   ''posterior_mode'', ')
            disp('                   ''posterior_mean'', ')
            disp('                   ''posterior_median'', ')
            disp('                   ''prior_mode'' or')
            disp('                   ''prior_mean''.')
            error('imcforecast:: Wrong argument type!')
        end
    else
        xparam = options_cond_fcst.parameter_set;
        if length(xparam)~=length(M_.params)
            error('imcforecast:: The dimension of the vector of parameters doesn''t match the number of estimated parameters!')
        end
    end
    set_parameters(xparam);
    [dataset_,dataset_info] = makedataset(options_);
    data = transpose(dataset_.data);
    data_index = dataset_info.missing.aindex;
    gend = dataset_.nobs;
    missing_value = dataset_info.missing.state;

    %store qz_criterium
    qz_criterium_old=options_.qz_criterium;
    options_=select_qz_criterium_value(options_);
    options_smoothed_state_uncertainty_old = options_.smoothed_state_uncertainty;
    [atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff,aK,T,R,P,PK,decomp,trend_addition,state_uncertainty,M_,oo_,options_,bayestopt_] = DsgeSmoother(xparam,gend,data,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
    options_.smoothed_state_uncertainty = options_smoothed_state_uncertainty_old;
    %get constant part
    if options_.noconstant
        constant = zeros(size(ys,1),options_cond_fcst.periods+1);
    else
        if options_.loglinear
            constant = repmat(log(ys),1,options_cond_fcst.periods+1);
        else
            constant = repmat(ys,1,options_cond_fcst.periods+1);
        end
    end
    %get trend part (which also takes care of prefiltering); needs to
    %include the last period
    if bayestopt_.with_trend == 1
        [trend_addition] =compute_trend_coefficients(M_,options_,size(bayestopt_.smoother_mf,1),gend+options_cond_fcst.periods);
        trend_addition = trend_addition(:,gend:end);
    else
        trend_addition=zeros(size(bayestopt_.smoother_mf,1),1+options_cond_fcst.periods);
    end
    % add trend to constant
    for obs_iter=1:length(options_.varobs)
        j = strmatch(options_.varobs{obs_iter},M_.endo_names,'exact');
        constant(j,:) = constant(j,:)+trend_addition(obs_iter,:);
    end
    trend = constant(oo_.dr.order_var,:);
    InitState(:,1) = atT(:,end);
else
    qz_criterium_old=options_.qz_criterium;
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    end
    graph_title='Calibration';
    if ~isfield(oo_.dr,'kstate')
        error('You need to call stoch_simul before conditional_forecast')
    end
end

if options_.logged_steady_state %if steady state was previously logged, undo this
    oo_.dr.ys=exp(oo_.dr.ys);
    oo_.steady_state=exp(oo_.steady_state);
    options_.logged_steady_state=0;
end

[T,R,ys,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);

if options_.loglinear && isfield(oo_.dr,'ys') && options_.logged_steady_state==0 %log steady state
    oo_.dr.ys=log_variable(1:M_.endo_nbr,oo_.dr.ys,M_);
    ys=oo_.dr.ys;
    oo_.steady_state=log_variable(1:M_.endo_nbr,oo_.steady_state,M_);
    options_.logged_steady_state=1; %set option for use in stoch_simul
end

if ~isdiagonal(M_.Sigma_e)
    warning(sprintf('The innovations are correlated (the covariance matrix has non zero off diagonal elements), the results of the conditional forecasts will\ndepend on the ordering of the innovations (as declared after varexo) because a Cholesky decomposition is used to factorize the covariance matrix.\n\n=> It is preferable to declare the correlations in the model block (explicitly imposing the identification restrictions), unless you are satisfied\nwith the implicit identification restrictions implied by the Cholesky decomposition.'))
    sQ = chol(M_.Sigma_e,'lower');
else
    sQ = sqrt(M_.Sigma_e);
end

if ~estimated_model
    if isempty(M_.endo_histval)
        y0 = ys;
    else
        if options_.loglinear
            %make sure that only states are updated (controls have value of 0 in vector)
            y0=zeros(size(ys));
            y0_logged = log_variable(1:M_.endo_nbr,M_.endo_histval,M_);
            y0(M_.endo_histval~=0)=y0_logged(M_.endo_histval~=0);
        else
            y0 = M_.endo_histval;
        end
    end
    InitState(:,1) = y0(oo_.dr.order_var)-ys(oo_.dr.order_var,:); %initial state in deviations from steady state
    trend = repmat(ys(oo_.dr.order_var,:),1,options_cond_fcst.periods+1); %trend needs to contain correct steady state
end




NumberOfStates = length(InitState);
FORCS1 = zeros(NumberOfStates,options_cond_fcst.periods+1,options_cond_fcst.replic);

FORCS1(:,1,:) = repmat(InitState,1,options_cond_fcst.replic); %set initial steady state to deviations from steady state in first period

EndoSize = M_.endo_nbr;
ExoSize = M_.exo_nbr;

n1 = size(constrained_vars,1);
n2 = size(options_cond_fcst.controlled_varexo,1);
constrained_vars(:,1)=oo_.dr.inv_order_var(constrained_vars); % must be in decision rule order

if n1 ~= n2
    error(['imcforecast:: The number of constrained variables doesn''t match the number of controlled shocks'])
end

idx = [];
jdx = [];

for i = 1:n1
    idx = [idx ; constrained_vars(i,:)];
    %     idx = [idx ; oo_.dr.inv_order_var(constrained_vars(i,:))];
    jdx = [jdx ; strmatch(deblank(options_cond_fcst.controlled_varexo(i,:)),M_.exo_names,'exact')];
end
mv = zeros(n1,NumberOfStates);
mu = zeros(ExoSize,n2);
for i=1:n1
    mv(i,idx(i)) = 1;
    mu(jdx(i),i) = 1;
end

% number of periods with constrained values
cL = size(constrained_paths,2);

%transform constrained periods into deviations from steady state; note that
%trend includes last actual data point and therefore we need to start in
%period 2
constrained_paths = bsxfun(@minus,constrained_paths,trend(idx,2:1+cL));

FORCS1_shocks = zeros(n1,cL,options_cond_fcst.replic);

%randn('state',0);

for b=1:options_cond_fcst.replic %conditional forecast using cL set to constrained values
    shocks = sQ*randn(ExoSize,options_cond_fcst.periods);
    shocks(jdx,:) = zeros(length(jdx),options_cond_fcst.periods);
    [FORCS1(:,:,b), FORCS1_shocks(:,:,b)] = mcforecast3(cL,options_cond_fcst.periods,constrained_paths,shocks,FORCS1(:,:,b),T,R,mv, mu);
    FORCS1(:,:,b)=FORCS1(:,:,b)+trend; %add trend
end

mFORCS1 = mean(FORCS1,3);
mFORCS1_shocks = mean(FORCS1_shocks,3);

tt = (1-options_cond_fcst.conditional_forecast.conf_sig)/2;
t1 = round(options_cond_fcst.replic*tt);
t2 = round(options_cond_fcst.replic*(1-tt));

forecasts.controlled_variables = constrained_vars;
forecasts.instruments = options_cond_fcst.controlled_varexo;

for i = 1:EndoSize
    forecasts.cond.Mean.(deblank(M_.endo_names(oo_.dr.order_var(i),:)))= mFORCS1(i,:)';
    tmp = sort(squeeze(FORCS1(i,:,:))');
    forecasts.cond.ci.(deblank(M_.endo_names(oo_.dr.order_var(i),:))) = [tmp(t1,:)' ,tmp(t2,:)' ]';
end

for i = 1:n1
    forecasts.controlled_exo_variables.Mean.(deblank(options_cond_fcst.controlled_varexo(i,:))) = mFORCS1_shocks(i,:)';
    tmp = sort(squeeze(FORCS1_shocks(i,:,:))');
    forecasts.controlled_exo_variables.ci.(deblank(options_cond_fcst.controlled_varexo(i,:))) = [tmp(t1,:)' ,tmp(t2,:)' ]';
end

clear FORCS1 mFORCS1_shocks;

FORCS2 = zeros(NumberOfStates,options_cond_fcst.periods+1,options_cond_fcst.replic);
FORCS2(:,1,:) = repmat(InitState,1,options_cond_fcst.replic); %set initial steady state to deviations from steady state in first period

%randn('state',0);

for b=1:options_cond_fcst.replic %conditional forecast using cL set to 0
    shocks = sQ*randn(ExoSize,options_cond_fcst.periods);
    shocks(jdx,:) = zeros(length(jdx),options_cond_fcst.periods);
    FORCS2(:,:,b) = mcforecast3(0,options_cond_fcst.periods,constrained_paths,shocks,FORCS2(:,:,b),T,R,mv, mu)+trend;
end

mFORCS2 = mean(FORCS2,3);

for i = 1:EndoSize
    forecasts.uncond.Mean.(deblank(M_.endo_names(oo_.dr.order_var(i),:)))= mFORCS2(i,:)';
    tmp = sort(squeeze(FORCS2(i,:,:))');
    forecasts.uncond.ci.(deblank(M_.endo_names(oo_.dr.order_var(i),:))) = [tmp(t1,:)' ,tmp(t2,:)' ]';
end
forecasts.graph.title=graph_title;
forecasts.graph.fname=M_.fname;

%reset qz_criterium
options_.qz_criterium=qz_criterium_old;

save('conditional_forecasts.mat','forecasts');