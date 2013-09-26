% computes the steady state of fs2000 analyticaly
% largely inspired by the program of F. Schorfheide

% Copyright (C) 2004-2010 Dynare Team
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

function [ys,check] = fs2000_steadystate(ys,exe)
  global M_
  
  alp = M_.params(1); 
  bet = M_.params(2); 
  gam = M_.params(3); 
  mst = M_.params(4); 
  rho = M_.params(5); 
  psi = M_.params(6); 
  del = M_.params(7); 

  check = 0;
  
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
  
  ys =[
m     
P     
c     
e     
W     
R     
k     
d     
n     
l     
gy_obs
gp_obs
y     
dA          ];
