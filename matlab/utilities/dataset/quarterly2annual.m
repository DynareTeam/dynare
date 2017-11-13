function [ya, yass, gya, gyass] = quarterly2annual(y,yss,GYTREND0,type,islog,aux)
% function [ya, yass, gya, gyass] = quarterly2annual(y,yss,GYTREND0,type,islog,aux)
% transforms quarterly (log-)level time series to annual level and growth rate
% it accounts for stock/flow/deflator series.
%
% INPUTS
% y        quarterly time series
% yss      steady state of y
% GYTREND0 growth rate of y
% type     1 sum (default)
%          2 average
%          3 last period (Q4)
%          4 geometric average
%          5 annual price as quantity weighted average
%          6 annual quantity from average price
%          7 annual nominal from Q real and deflator
% islog    0 level (default)
%          1 log-level
%          2 growth rate Q frequency
% aux      optional input used when type>4
%
%
% OUTPUTS
% ya       annual (log-)level
% yass     annual steadystate (log-)level
% gya      annual growth rate
% gyass    annual growth rate steadystate

% Copyright (C) 2017 Dynare Team
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

if nargin ==0
    disp('[ya, yass, gya, gyass] = quarterly2annual(y,yss,GYTREND0,type,islog);')
    return
end

if nargin<4 || isempty(type)
    type=1;
end
if nargin<5 || isempty(islog)
    islog=0;
end
if isstruct(aux)
    yaux=aux.y;
    yauxss=aux.yss;
    islogaux=aux.islog;
    GYTREND0aux=aux.GYTREND0;
    typeaux=aux.type;
    if islogaux
        yaux=exp(yaux+yauxss);
        yauxss=exp(yauxss);
        yaux=yaux-yauxss;
    end
elseif type > 4
    error('TYPE>4 requires auxiliary variable!')
end
if islog == 2
    % construct loglevel out of growth rate
    y = cumsum(y);
    yss=0;
    islog=1;
end
if islog == 1
    y=exp(y+yss);
    yss=exp(yss);
    y=y-yss;
end
switch type
  case 1
    yass = yss*(exp(-GYTREND0*3)+exp(-GYTREND0*2)+exp(-GYTREND0)+1);
    tmp = lagged(y,3)*exp(-GYTREND0*3)+lagged(y,2)*exp(-GYTREND0*2)+lagged(y,1)*exp(-GYTREND0)+y; % annualized level
    ya = tmp(4:4:end);
  case 2
    yass = yss*(exp(-GYTREND0*3)+exp(-GYTREND0*2)+exp(-GYTREND0)+1)/4;
    tmp = (lagged(y,3)*exp(-GYTREND0*3)+lagged(y,2)*exp(-GYTREND0*2)+lagged(y,1)*exp(-GYTREND0)+y)/4; % annualized level
    ya = tmp(4:4:end);
  case 3
    yass=yss;
    tmp = y;
    ya = tmp(4:4:end);
  case 4
    yass = yss*(exp(-GYTREND0*3/2));
    tmp = (lagged(y+yss,3)*exp(-GYTREND0*3).*lagged(y+yss,2)*exp(-GYTREND0*2).*lagged(y+yss,1)*exp(-GYTREND0).*(y+yss)).^(1/4); % annualized level
    tmp = tmp - yass;
    ya = tmp(4:4:end);
  case 5
    % nominal series
    yn = (y+yss).*(yaux+yauxss) - yss.*yauxss;
    [yna, ynass] = quarterly2annual(yn,yss.*yauxss,GYTREND0+GYTREND0aux,typeaux,0,0);
    % real series
    [yra, yrass] = quarterly2annual(yaux,yauxss,GYTREND0aux,typeaux,0,0);
    % deflator
    yass = ynass/yrass;
    ya = (yna+ynass)./(yra+yrass)-yass;
  case 6
    % nominal series
    yn = (y+yss).*(yaux+yauxss) - yss.*yauxss;
    [yna, ynass] = quarterly2annual(yn,yss.*yauxss,GYTREND0+GYTREND0aux,typeaux,0,0);
    % deflator
    [pa, pass] = quarterly2annual(yaux,yauxss,GYTREND0aux,2,0,0);
    % real series
    yass = ynass/pass;
    ya = (yna+ynass)./(pa+pass)-yass;
  case 7
    % nominal series
    yn = (y+yss).*(yaux+yauxss) - yss.*yauxss;
    [ya, yass] = quarterly2annual(yn,yss.*yauxss,GYTREND0+GYTREND0aux,typeaux,0,0);
    GYTREND0=GYTREND0+GYTREND0aux;
  otherwise
    error('Wrong type input')
end

% annual growth rate
gyass = GYTREND0*4;
gya = (ya+yass)./(lagged(ya,1)+yass).*exp(4*GYTREND0)-1-gyass;

if islog
    ya=log(ya+yass);
    yass=log(yass);
    ya=ya-yass;
end
