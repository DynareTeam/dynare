function imcforecast(constrained_paths, constrained_vars, options_cond_fcst, constrained_perfect_foresight)
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

% Copyright (C) 2006-2013 Dynare Team
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

global options_ oo_ M_ bayestopt_


if isfield(options_cond_fcst, 'simulation_type')
    if strcmp(options_cond_fcst.simulation_type, 'deterministic')
        disp('deterministic condtional forecast');
        det_cond_forecast(constrained_paths, constrained_vars, options_cond_fcst, constrained_perfect_foresight);
        return;
    end
end
    
if ~isfield(options_cond_fcst,'parameter_set') || isempty(options_cond_fcst.parameter_set)
    options_cond_fcst.parameter_set = 'posterior_mode';
end

if ~isfield(options_cond_fcst,'replic') || isempty(options_cond_fcst.replic)
    options_cond_fcst.replic = 5000;
end

if ~isfield(options_cond_fcst,'periods') || isempty(options_cond_fcst.periods)
    options_cond_fcst.periods = 40;
end

if ~isfield(options_cond_fcst,'conf_sig') || isempty(options_cond_fcst.conf_sig)
    options_cond_fcst.conf_sig = .8;
end

if isequal(options_cond_fcst.parameter_set,'calibration')
    estimated_model = 0;
else
    estimated_model = 1;
end

if estimated_model
    if ischar(options_cond_fcst.parameter_set)
        switch options_cond_fcst.parameter_set
          case 'posterior_mode'
            xparam = get_posterior_parameters('mode');
            graph_title='Posterior Mode';
          case 'posterior_mean'
            xparam = get_posterior_parameters('mean');
            graph_title='Posterior Mean';
          case 'posterior_median'
            xparam = get_posterior_parameters('median');
            graph_title='Posterior Median';
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

    % Load and transform data.
    transformation = [];
    if options_.loglinear && ~options_.logdata
        transformation = @log;
    end
    xls.sheet = options_.xls_sheet;
    xls.range = options_.xls_range;
    
    if ~isfield(options_,'nobs')
        options_.nobs = [];
    end
    
    dataset_ = initialize_dataset(options_.datafile,options_.varobs,options_.first_obs,options_.nobs,transformation,options_.prefilter,xls);

    data = dataset_.data;
    data_index = dataset_.missing.aindex;
    gend = options_.nobs;
    missing_value = dataset_.missing.state;

    [atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(xparam,gend,data,data_index,missing_value);

    trend = repmat(ys,1,options_cond_fcst.periods+1);
    for i=1:M_.endo_nbr
        j = strmatch(deblank(M_.endo_names(i,:)),options_.varobs,'exact');
        if ~isempty(j)
            trend(i,:) = trend(i,:)+trend_coeff(j)*(gend+(0:options_cond_fcst.periods));
        end
    end
    trend = trend(oo_.dr.order_var,:);

    InitState(:,1) = atT(:,end);
else
    InitState(:,1) = zeros(M_.endo_nbr,1);
    trend = repmat(oo_.steady_state(oo_.dr.order_var),1,options_cond_fcst.periods+1);
    graph_title='Calibration';
end

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end
[T,R,ys,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);

sQ = sqrt(M_.Sigma_e);

NumberOfStates = length(InitState);
FORCS1 = zeros(NumberOfStates,options_cond_fcst.periods+1,options_cond_fcst.replic);

FORCS1(:,1,:) = repmat(InitState,1,options_cond_fcst.replic);

EndoSize = M_.endo_nbr;
ExoSize = M_.exo_nbr;

n1 = size(constrained_vars,1);
n2 = size(options_cond_fcst.controlled_varexo,1);

if n1 ~= n2
    error(['imcforecast:: The number of constrained variables doesn''t match the number of controlled shocks'])
end

idx = [];
jdx = [];

for i = 1:n1
    idx = [idx ; oo_.dr.inv_order_var(constrained_vars(i,:))];
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

constrained_paths = bsxfun(@minus,constrained_paths,trend(idx,1:cL));

%randn('state',0);

for b=1:options_cond_fcst.replic
    shocks = sQ*randn(ExoSize,options_cond_fcst.periods);
    shocks(jdx,:) = zeros(length(jdx),options_cond_fcst.periods);
    FORCS1(:,:,b) = mcforecast3(cL,options_cond_fcst.periods,constrained_paths,shocks,FORCS1(:,:,b),T,R,mv, mu)+trend;
end

mFORCS1 = mean(FORCS1,3);

tt = (1-options_cond_fcst.conf_sig)/2;
t1 = round(options_cond_fcst.replic*tt);
t2 = round(options_cond_fcst.replic*(1-tt));

forecasts.controled_variables = constrained_vars;
forecasts.instruments = options_cond_fcst.controlled_varexo;

for i = 1:EndoSize
    eval(['forecasts.cond.Mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = mFORCS1(i,:)'';']);
    tmp = sort(squeeze(FORCS1(i,:,:))');
    eval(['forecasts.cond.ci.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ...
          ' = [tmp(t1,:)'' ,tmp(t2,:)'' ]'';']);
end

clear FORCS1;

FORCS2 = zeros(NumberOfStates,options_cond_fcst.periods+1,options_cond_fcst.replic);
for b=1:options_cond_fcst.replic
    FORCS2(:,1,b) = InitState;
end

%randn('state',0);

for b=1:options_cond_fcst.replic
    shocks = sQ*randn(ExoSize,options_cond_fcst.periods);
    shocks(jdx,:) = zeros(length(jdx),options_cond_fcst.periods);
    FORCS2(:,:,b) = mcforecast3(0,options_cond_fcst.periods,constrained_paths,shocks,FORCS2(:,:,b),T,R,mv, mu)+trend;
end

mFORCS2 = mean(FORCS2,3);

for i = 1:EndoSize
    eval(['forecasts.uncond.Mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = mFORCS2(i,:)'';']);
    tmp = sort(squeeze(FORCS2(i,:,:))');
    eval(['forecasts.uncond.ci.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ...
          ' = [tmp(t1,:)'' ,tmp(t2,:)'' ]'';']);
end
forecasts.graph.title=graph_title;
forecasts.graph.fname=M_.fname;

save('conditional_forecasts.mat','forecasts');