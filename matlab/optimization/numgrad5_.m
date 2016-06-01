function [g, badg] = numgrad5_(fcn,f0,x,penalty,epsilon,varargin)
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

rescale_step_length = 0;

delta = epsilon;
n = length(x);
g = zeros(n,1);

badg = 0;

if rescale_step_length
    scale = [];
else
    scale = ones(n,1);
end

for i=1:n
    xiold = x(i);
    h = step_length_correction(xiold,scale,i)*delta;
    x(i) = xiold+h;
    [f1,junk1,cost_flag1,] = penalty_objective_function(x, fcn, penalty, varargin{:});
    x(i) = xiold-h;
    [f2,junk1,cost_flag2] = penalty_objective_function(x, fcn, penalty, varargin{:});
    x(i) = xiold+2*h;
    [f3,junk1,cost_flag3] = penalty_objective_function(x, fcn, penalty, varargin{:});
    x(i) = xiold-2*h;
    [f4,junk1,cost_flag4] = penalty_objective_function(x, fcn, penalty, varargin{:});
    if f0<f1 && f1<f3 && f0<f2 && f2<f4
        g0 = 0;
    else
        g0 = (8*(f1 - f2)+ f4-f3) / (12*h);
    end
    if abs(g0)< 1e15
        g(i) = g0;
    else
        disp('Bad gradient -----------------------------------')
        fprintf('Gradient w.r.t. parameter number %3d (x=%16.8f,h=%16.8f,g0=%16.8f)\n',i,xiold,h,g0)
        g(i) = 0;
        badg = 1;
    end
end