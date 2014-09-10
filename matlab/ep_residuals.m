function r = ep_residuals(x, y, ix, iy, steadystate, dr, maximum_lag, endo_nbr)
% Inversion of the extended path simulation approach. This routine computes the innovations needed to
% reproduce the time path of a subset of endogenous variables.    
%    
% INPUTS
%  o x    [double]   n*1 vector, time t innovations.    
%  o y    [double]   n*1 vector, time t restricted endogenous variables.
%  o ix   [integer]  index of control innovations in the full vector of innovations.
%  o iy   [integer]  index of controlled variables in the full vector of endogenous variables.    
%  o s    [double]   m*1 vector, endogenous variables at time t-1. 
%    
% 
% OUTPUTS
%  o r    [double]  n*1 vector of residuals.
%    
% ALGORITHM
%  
% SPECIAL REQUIREMENTS

% Copyright (C) 2010 Dynare Team.
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

global oo_ options_

persistent k1 k2 weight

if isempty(k1)
    k1 = [maximum_lag:-1:1];
    k2 = dr.kstate(find(dr.kstate(:,2) <= maximum_lag+1),[1 2]);
    k2 = k2(:,1)+(maximum_lag+1-k2(:,2))*endo_nbr;
    weight = 0.0;
end

verbose = options_.ep.verbosity;

% Copy the shocks in exo_simul.
oo_.exo_simul(maximum_lag+1,ix) = exp(transpose(x));
exo_simul = log(oo_.exo_simul);

% Compute the initial solution path for the endogenous variables using a first order approximation.
if verbose
    disp('ep_residuals:: Set initial condition for endogenous variable paths.')
end
initial_path = oo_.endo_simul;
for i = maximum_lag+1:size(oo_.exo_simul)
    tempx1 = oo_.endo_simul(dr.order_var,k1);
    tempx2 = bsxfun(@minus,tempx1,dr.ys(dr.order_var));
    tempx = tempx2(k2);
    initial_path(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx2(k2)+dr.ghu*transpose(exo_simul(i,:));
    k1 = k1+1;
end
oo_.endo_simul = weight*initial_path + (1-weight)*oo_.endo_simul;

info = perfect_foresight_simulation(dr,steadystate);
if verbose>1
    info
    info.iterations.errors
end

r = y-transpose(oo_.endo_simul(maximum_lag+1,iy));

%(re)Set k1 (indices for the initial conditions)
k1 = [maximum_lag:-1:1];