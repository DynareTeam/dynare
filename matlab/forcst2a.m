function yf=forcst2a(y0,dr,e)
% function yf=forcst2a(y0,dr,e)
% computes forecasts based on first order model solution, assuming the absence of shocks 
% Inputs:
%   - y0        [endo_nbr by maximum_endo_lag]          matrix of starting values
%   - dr        [structure]                             structure with Dynare decision rules
%   - e         [horizon by exo_nbr]                    matrix with shocks
%
% Outputs:
%   - yf        [horizon+maximum_endo_lag,endo_nbr]               matrix of forecasts
%
% Copyright (C) 2008-2015 Dynare Team
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

global M_ options_

endo_nbr = M_.endo_nbr;
ykmin_ = M_.maximum_endo_lag;

horizon = size(e,1);

k1 = [ykmin_:-1:1];
k2 = dr.kstate(find(dr.kstate(:,2) <= ykmin_+1),[1 2]);
k2 = k2(:,1)+(ykmin_+1-k2(:,2))*endo_nbr;

yf = zeros(horizon+ykmin_,endo_nbr);
yf(1:ykmin_,:) = y0';

j = ykmin_*endo_nbr;
for i=ykmin_+(1:horizon)
    tempx = yf(k1,:)';
    yf(i,:) = tempx(k2)'*dr.ghx';
    k1 = k1+1;
end

yf(:,dr.order_var) = yf;

