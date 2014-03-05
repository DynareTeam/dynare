function forcs = mcforecast3(cL,H,mcValue,shocks,forcs,T,R,mv,mu)
% forcs = mcforecast3(cL,H,mcValue,shocks,forcs,T,R,mv,mu)
% Computes the shock values for constrained forecasts necessary to keep
% endogenous variables at their constrained paths
% 
% INPUTS
%  o cL             [scalar]                            number of controlled periods
%  o H              [scalar]                            number of forecast periods
%  o mcValue        [n_controlled_vars by cL double]    paths for constrained variables
%  o shocks         [nexo by H double]      shock values draws (with zeros for controlled_varexo)
%  o forcs
%  o T              [n_endovars by n_endovars double]       transition matrix of the state equation.
%  o R              [n_endovars by n_exo double]           matrix relating the endogenous variables to the innovations in the state equation.
%  o mv             [n_controlled_exo by n_endovars boolean]        indicator vector  selecting constrained endogenous variables
%  o mu             [n_controlled_vars by nexo boolean]             indicator vector
%                                                   selecting controlled exogenous variables
% 
% Algorithm:
%   Relies on state-space form:
%       y_t=T*y_{t-1}+R*shocks(:,t)
%   Shocks are split up into shocks_uncontrolled and shockscontrolled while
%   the endogenous variables are also split up into controlled and
%   uncontrolled ones to get:
%       y_t(controlled_vars_index)=T*y_{t-1}(controlled_vars_index)+R(controlled_vars_index,uncontrolled_shocks_index)*shocks_uncontrolled_t
%                    + R(controlled_vars_index,controlled_shocks_index)*shocks_controlled_t
%
%   This is then solved to get:
%       shocks_controlled_t=(y_t(controlled_vars_index)-(T*y_{t-1}(controlled_vars_index)+R(controlled_vars_index,uncontrolled_shocks_index)*shocks_uncontrolled_t)/R(controlled_vars_index,controlled_shocks_index)
%
%   After obtaining the shocks, and for uncontrolled periods, the state-space representation
%       y_t=T*y_{t-1}+R*shocks(:,t)
%   is used for forecasting
%  
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

if cL
    e = zeros(size(mcValue,1),cL);
    for t = 1:cL
        e(:,t) = inv(mv*R*mu)*(mcValue(:,t)-mv*T*forcs(:,t)-mv*R*shocks(:,t));
        forcs(:,t+1) = T*forcs(:,t)+R*(mu*e(:,t)+shocks(:,t));
    end
end
for t = cL+1:H
    forcs(:,t+1) = T*forcs(:,t)+R*shocks(:,t);
end