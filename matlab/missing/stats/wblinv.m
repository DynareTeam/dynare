function t = wblinv(proba, scale, shape)   % --*-- Unitary tests --*--

% Inverse cumulative distribution function.
%
% INPUTS 
% - proba [double] Probability, scalar between 0 and 1.
% - scale [double] Positive hyperparameter.
% - shape [double] Positive hyperparameter.
%
% OUTPUTS 
% - t     [double] scalar such that P(X<=t)=proba

% Copyright (C) 2015 Dynare Team
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

% Check input arguments.

if nargin<3
    error('Three input arguments required!')
end

if ~isnumeric(proba) || ~isscalar(proba) || ~isreal(proba) || proba<0 || proba>1
    error('First input argument must be a real scalar between 0 and 1 (probability)!')
end

if ~isnumeric(scale) || ~isscalar(scale) || ~isreal(scale) || scale<=0
    error('Second input argument must be a real positive scalar (scale parameter of the Weibull distribution)!')
end

if ~isnumeric(shape) || ~isscalar(shape) || ~isreal(shape) || shape<=0
    error('Third input argument must be a real positive scalar (shape parameter of the Weibull distribution)!')
end


if proba<2*eps()
    t = 0;
    return
end

if proba>1-2*eps()
    t = Inf;
    return
end

t = exp(log(scale)+log(-log(1-proba))/shape);

%@test:1
%$ try
%$    x = wblinv(0, 1, 2);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = isequal(x, 0);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$    x = wblinv(1, 1, 2);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = isinf(x);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ scales = [.5, 1, 5];
%$ shapes = [.1, 1, 2];
%$ x = NaN(9,1);
%$
%$ try
%$    k = 0;
%$    for i=1:3
%$       for j=1:3
%$           k = k+1;
%$           x(k) = wblinv(.5, scales(i), shapes(j));
%$       end
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    k = 1;
%$    for i=1:3
%$       for j=1:3
%$           k = k+1;
%$           t(k) = abs(x(k-1)-scales(i)*log(2)^(1/shapes(j)))<1e-12;
%$       end
%$    end
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ debug = false;
%$ scales = [ .5, 1, 5];
%$ shapes = [ 1, 2, 3];
%$ x = NaN(9,1);
%$ p = 1e-1;
%$
%$ try
%$    k = 0;
%$    for i=1:3
%$       for j=1:3
%$           k = k+1;
%$           x(k) = wblinv(p, scales(i), shapes(j));
%$       end
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    k = 1;
%$    for i=1:3
%$       for j=1:3
%$           k = k+1;
%$           shape = shapes(j);
%$           scale = scales(i);
%$           density = @(z) exp(lpdfgweibull(z,shape,scale));
%$           if debug
%$               [shape, scale, x(k-1)]
%$           end
%$           if isoctave
%$               s = quadv(density, 0, x(k-1),1e-10);
%$           else
%$               s = integral(density, 0, x(k-1));
%$           end
%$           if debug
%$               [s, abs(p-s)]
%$           end
%$         if isoctave
%$           t(k) = abs(p-s)<1e-10;  
%$         else
%$           t(k) = abs(p-s)<1e-12;
%$         end
%$       end
%$    end
%$ end
%$ T = all(t);
%@eof:4

