function [forecast,info] = dyn_forecast(var_list,M,options,oo,task,dataset_info)
% function dyn_forecast(var_list,M,options,oo,task,dataset_info)
%   computes mean forecast for a given value of the parameters
%   compues also confidence band for the forecast    
%
% INPUTS
%   var_list:    list of variables (character matrix)
%   M:           Dynare model structure
%   options:     Dynare options structure
%   oo:          Dynare results structure
%   task:        indicates how to initialize the forecast
%                either 'simul' or 'smoother'
%   dataset_info:   Various informations about the dataset (descriptive statistics and missing observations).

% OUTPUTS
%   nothing is returned but the procedure saves output
%   in oo_.forecast.Mean
%      oo_.forecast.HPDinf
%      oo_.forecast.HPDsup
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2016 Dynare Team
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

if nargin<6 && options.prefilter
    error('The prefiltering option is not allowed without providing a dataset')
elseif nargin==6
    mean_varobs=dataset_info.descriptive.mean';
end
    
info = 0;

oo=make_ex_(M,options,oo);

maximum_lag = M.maximum_lag;

endo_names = M.endo_names;
if isempty(var_list)
    var_list = endo_names(1:M.orig_endo_nbr, :);
end
i_var = [];
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t an endogenous variable'])
    end
    i_var = [i_var; tmp];
end

n_var = length(i_var);

trend = 0;
switch task
  case 'simul'
    horizon = options.periods;
    if horizon == 0
        horizon = 5;
    end
    if isempty(M.endo_histval)
        if options.loglinear && ~options.logged_steady_state
            y0 = repmat(log(oo.dr.ys),1,maximum_lag);
        else
            y0 = repmat(oo.dr.ys,1,maximum_lag);
        end
    else
        if options.loglinear
            y0 = log_variable(1:M.endo_nbr,M.endo_histval,M);
        else
            y0 = M.endo_histval;
        end
    end
  case 'smoother'
    horizon = options.forecast;
    y_smoothed = oo.SmoothedVariables;
    y0 = zeros(M.endo_nbr,maximum_lag);
    for i = 1:M.endo_nbr
        v_name = deblank(M.endo_names(i,:));
        y0(i,:) = y_smoothed.(v_name)(end-maximum_lag+1:end); %includes steady state or mean, but simult_ will subtract only steady state
        % 2. Subtract mean/steady state and add steady state; takes care of prefiltering        
        if isfield(oo.Smoother,'Constant') && isfield(oo.Smoother.Constant,v_name)
            y0(i,:)=y0(i,:)-oo.Smoother.Constant.(v_name)(end-maximum_lag+1:end); %subtract mean or steady state
            if options.loglinear
                y0(i,:)=y0(i,:)+log_variable(i,oo.dr.ys,M);
            else
                y0(i,:)=y0(i,:)+oo.dr.ys(strmatch(v_name,deblank(M.endo_names),'exact'));
            end
        end                
        % 2. Subtract trend
        if isfield(oo.Smoother,'Trend') && isfield(oo.Smoother.Trend,v_name)
            y0(i,:)=y0(i,:)-oo.Smoother.Trend.(v_name)(end-maximum_lag+1:end); %subtract trend, which is not subtracted by simult_
        end
    end
    gend = options.nobs;
    if isfield(oo.Smoother,'TrendCoeffs')
        var_obs = options.varobs;
        endo_names = M.endo_names;
        order_var = oo.dr.order_var;
        i_var_obs = [];
        trend_coeffs = [];
        for i=1:length(var_obs)
            tmp = strmatch(var_obs{i},endo_names(i_var,:),'exact');
            trend_var_index=strmatch(var_obs{i},M.endo_names,'exact');
            if ~isempty(tmp)
                i_var_obs = [ i_var_obs; tmp];
                trend_coeffs = [trend_coeffs; oo.Smoother.TrendCoeffs(trend_var_index)];
            end
        end
        if ~isempty(trend_coeffs) 
            trend = trend_coeffs*(options.first_obs+gend-1+(1-M.maximum_lag:horizon)); 
            if options.prefilter
                trend = trend - repmat(mean(trend_coeffs*[options.first_obs:options.first_obs+gend-1],2),1,horizon+1); %subtract mean trend
            end
        end
    else
        trend_coeffs=zeros(length(options.varobs),1);
    end
  otherwise
    error('Wrong flag value')
end 

if M.exo_det_nbr == 0
    if isequal(M.H,0)
        [yf,int_width] = forcst(oo.dr,y0,horizon,var_list,M,oo,options);
    else
        [yf,int_width,int_width_ME] = forcst(oo.dr,y0,horizon,var_list,M,oo,options);        
    end
else
    exo_det_length = size(oo.exo_det_simul,1)-M.maximum_lag;
    if horizon > exo_det_length
        ex = zeros(horizon,M.exo_nbr);
        oo.exo_det_simul = [ oo.exo_det_simul;...
                            repmat(oo.exo_det_steady_state',...
                                   horizon- ... 
                                   exo_det_length,1)];
    elseif horizon <= exo_det_length 
        ex = zeros(exo_det_length,M.exo_nbr); 
    end
    if isequal(M.H,0)
        [yf,int_width] = simultxdet(y0,ex,oo.exo_det_simul,...
                                options.order,var_list,M,oo,options);
    else
        [yf,int_width,int_width_ME] = simultxdet(y0,ex,oo.exo_det_simul,...
                                options.order,var_list,M,oo,options);
    end
end

if ~isscalar(trend) %add trend back to forecast
    yf(i_var_obs,:) = yf(i_var_obs,:) + trend;
end

if options.loglinear == 1
    if options.prefilter == 1 %subtract steady state and add mean for observables
        yf(i_var_obs,:)=yf(i_var_obs,:)-repmat(log(oo.dr.ys(i_var_obs)),1,horizon+M.maximum_lag)+ repmat(mean_varobs,1,horizon+M.maximum_lag);
    end
else
    if options.prefilter == 1 %subtract steady state and add mean for observables
        yf(i_var_obs,:)=yf(i_var_obs,:)-repmat(oo.dr.ys(i_var_obs),1,horizon+M.maximum_lag)+ repmat(mean_varobs,1,horizon+M.maximum_lag);
    end    
end

for i=1:n_var
    vname = deblank(var_list(i,:));
    forecast.Mean.(vname) = yf(i,maximum_lag+(1:horizon))';
    forecast.HPDinf.(vname)= yf(i,maximum_lag+(1:horizon))' - int_width(1:horizon,i);
    forecast.HPDsup.(vname) = yf(i,maximum_lag+(1:horizon))' + int_width(1:horizon,i);
    if ~isequal(M.H,0) && ismember(var_list(i,:),options.varobs)
        forecast.HPDinf_ME.(vname)= yf(i,maximum_lag+(1:horizon))' - int_width_ME(1:horizon,i);
        forecast.HPDsup_ME.(vname) = yf(i,maximum_lag+(1:horizon))' + int_width_ME(1:horizon,i);
    end
end

for i=1:M.exo_det_nbr
    forecast.Exogenous.(deblank(M.exo_det_names(i,:))) = oo.exo_det_simul(maximum_lag+(1:horizon),i);
end

if options.nograph == 0
    oo.forecast = forecast;
    forecast_graphs(var_list,M, oo,options)
end
