function [g, badg] = numgrad5(fcn,f0,x,penalty,epsilon,varargin)
% Computes the gradient of the objective function fcn using a five points
% formula if possible.
%
% Adapted from Sims' numgrad.m routine.
%
% See section 25.3.6 Abramovitz and Stegun (1972, Tenth Printing, December) Handbook of Mathematical Functions.
% http://www.math.sfu.ca/~cbm/aands/
%
% TODO Try Four points formula when cost_flag3=0 or cost_flag4=0.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/numgrad.m

% Copyright (C) 1993-2007 Christopher Sims
% Copyright (C) 2008-2016 Dynare Team

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

h = epsilon;
n = length(x);
g = zeros(n,1);
H = 12*h;

badg = 0;

for i=1:n
    xiold = x(i);
    x(i) = xiold+h;
    f1 = penalty_objective_function(x, fcn, penalty, varargin{:});
    x(i) = xiold-h;
    f2 = penalty_objective_function(x, fcn, penalty, varargin{:});
    x(i) = xiold+2*h;
    f3 = penalty_objective_function(x, fcn, penalty, varargin{:});
    x(i) = xiold-2*h;
    f4 = penalty_objective_function(x, fcn, penalty, varargin{:});
    g0 = (8*(f1-f2)+f4-f3)/H;
    if abs(g0)< 1e15
        g(i) = g0;
    else
        g(i) = 0;
        badg = 1;
    end
    x(i) = xiold;
end