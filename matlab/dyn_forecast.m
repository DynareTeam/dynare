function info = dyn_forecast(var_list,task,dataset_info)
% function dyn_forecast(var_list,task,dataset_info))
%   computes mean forecast for a given value of the parameters
%   computes also confidence band for the forecast    
%
% INPUTS
%   var_list:       list of variables (character matrix)
%   task:           indicates how to initialize the forecast
%                   either 'simul' or 'smoother'
%   dataset_info:   Various informations about the dataset (descriptive statistics and missing observations).
% OUTPUTS
%   nothing is returned but the procedure saves output
%   in oo_.forecast.Mean
%      oo_.forecast.HPDinf
%      oo_.forecast.HPDsup
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2015 Dynare Team
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

global options_ oo_ M_

if nargin<3 && options_.prefilter
    error('The prefiltering option is not allow without providing a dataset')
else
    mean_varobs=dataset_info.descriptive.mean';
end
    
info = 0;

make_ex_;

maximum_lag = M_.maximum_lag;

endo_names = M_.endo_names;
if isempty(var_list)
    var_list = endo_names(1:M_.orig_endo_nbr, :);
end
i_var = [];
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t and endogenous variable'])
    end
    i_var = [i_var; tmp];
end

n_var = length(i_var);

trend = 0;
switch task
  case 'simul'
    horizon = options_.periods;
    if horizon == 0
        horizon = 5;
    end
    if isempty(M_.endo_histval)
        y0 = repmat(oo_.dr.ys,1,maximum_lag);
    else
        y0 = M_.endo_histval;
    end
  case 'smoother'
    horizon = options_.forecast;
    y_smoothed = oo_.SmoothedVariables;
    y0 = zeros(M_.endo_nbr,maximum_lag);
    for i = 1:M_.endo_nbr
        v_name = deblank(M_.endo_names(i,:));
        y0(i,:) = y_smoothed.(v_name)(end-maximum_lag+1:end); %includes steady state or mean, but simult_ will subtract only steady state
        % 2. Subtract mean/steady state and add steady state; takes care of prefiltering        
        if isfield(oo_.Smoother,'Constant') && isfield(oo_.Smoother.Constant,v_name)
            y0(i,:)=y0(i,:)-oo_.Smoother.Constant.(v_name)(end-maximum_lag+1:end); %subtract mean or steady state
            if options_.loglinear
                y0(i,:)=y0(i,:)+log(oo_.dr.ys(strmatch(v_name,deblank(M_.endo_names),'exact')));
            else
                y0(i,:)=y0(i,:)+oo_.dr.ys(strmatch(v_name,deblank(M_.endo_names),'exact'));                
            end
        end                
        % 2. Subtract trend
        if isfield(oo_.Smoother,'Trend') && isfield(oo_.Smoother.Trend,v_name)
            y0(i,:)=y0(i,:)-oo_.Smoother.Trend.(v_name)(end-maximum_lag+1:end); %subtract trend, which is not subtracted by simult_
        end
    end
    gend = options_.nobs;
    if isfield(oo_.Smoother,'TrendCoeffs')
        var_obs = options_.varobs;
        endo_names = M_.endo_names;
        order_var = oo_.dr.order_var;
        i_var_obs = [];
        trend_coeffs = [];
        for i=1:length(var_obs)
            tmp = strmatch(var_obs{i},endo_names(i_var,:),'exact');
            if ~isempty(tmp)
                i_var_obs = [ i_var_obs; tmp];
                trend_coeffs = [trend_coeffs; oo_.Smoother.TrendCoeffs(i)];
            end
        end
        if ~isempty(trend_coeffs) 
            trend = trend_coeffs*(options_.first_obs+gend-1+(1-M_.maximum_lag:horizon)); 
            if options_.prefilter
                trend = trend - repmat(mean(trend_coeffs*[options_.first_obs:options_.first_obs+gend-1],2),1,horizon+1); %subtract mean trend
            end
        end
    else
        trend_coeffs=zeros(length(options_.varobs),1);
    end
  otherwise
    error('Wrong flag value')
end 

if M_.exo_det_nbr == 0
    [yf,int_width] = forcst(oo_.dr,y0,horizon,var_list);
else
    exo_det_length = size(oo_.exo_det_simul,1)-M_.maximum_lag;
    if horizon > exo_det_length
        ex = zeros(horizon,M_.exo_nbr);
        oo_.exo_det_simul = [ oo_.exo_det_simul;...
                            repmat(oo_.exo_det_steady_state',...
                                   horizon- ... 
                                   exo_det_length,1)];
    elseif horizon < exo_det_length 
        ex = zeros(exo_det_length,M_.exo_nbr); 
    end
    [yf,int_width] = simultxdet(y0,ex,oo_.exo_det_simul,...
                                options_.order,var_list,M_,oo_,options_);
end

if ~isscalar(trend) %add trend back to forecast
    yf(i_var_obs,:) = yf(i_var_obs,:) + trend;
end

if options_.loglinear == 1
    if options_.prefilter == 1 %subtract steady state and add mean for observables
        yf(i_var_obs,:)=yf(i_var_obs,:)-repmat(log(oo_.dr.ys(i_var_obs)),1,horizon+M_.maximum_lag)+ repmat(mean_varobs,1,horizon+M_.maximum_lag);
    end
else
    if options_.prefilter == 1 %subtract steady state and add mean for observables
        yf(i_var_obs,:)=yf(i_var_obs,:)-repmat(oo_.dr.ys(i_var_obs),1,horizon+M_.maximum_lag)+ repmat(mean_varobs,1,horizon+M_.maximum_lag);
    end    
end

for i=1:n_var
    eval(['oo_.forecast.Mean.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+(1:horizon))'';']);
    eval(['oo_.forecast.HPDinf.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+(1:horizon))''-' ...
          ' int_width(1:horizon,' int2str(i) ');']);
    eval(['oo_.forecast.HPDsup.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+(1:horizon))''+' ...
          ' int_width(1:horizon,' int2str(i) ');']);
end

for i=1:M_.exo_det_nbr
    eval(['oo_.forecast.Exogenous.' M_.exo_det_names(i,:) '= oo_.exo_det_simul(maximum_lag+(1:horizon),' int2str(i) ');']);
end

if options_.nograph == 0
    forecast_graphs(var_list,M_, oo_,options_)
end
