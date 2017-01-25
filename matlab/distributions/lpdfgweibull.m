function [ldens,Dldens,D2ldens] = lpdfgweibull(x,a,b,c)  % --*-- Unitary tests --*--

% Evaluates the logged Weibull PDF at x. 
%
% INPUTS 
% - x       [double]  m*n matrix of points where the (logged) density will be evaluated,
% - a       [double]  m*n matrix of First Weibull distribution parameters (shape parameter, k),
% - b       [double]  m*n matrix of Second Weibull distribution parameters (scale parameter, λ),
% - c       [double]  m*n matrix of Third Weibull distribution parameters (location parameter, default is 0).
%
% OUTPUTS 
% - ldens   [double]  m*n matrix of logged (generalized) Weibull densities.
% - Dldens  [double]  m*n matrix (first order derivatives w.r.t. x)
% - D2ldens [double]  m*n matrix (second order derivatives w.r.t. x)
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2015-2016 Dynare Team
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

% Initialize output arguments
ldens   = -Inf(size(x));
if nargout>1
    Dldens  = NaN(size(x));
    D2ldens = NaN(size(x));
end

% Check the number of input arguments
if nargin<3
    error('CheckInputArgs','At least three input arguments required!')
end

% Set default value for location parameter(s).
if nargin<4
    c = zeros(size(x));
end

% Reshape inputs if needed (and possible)
if ~isscalar(x)
    if isscalar(a)
        a = repmat(a, size(x));
    end 
    if isscalar(b)
        b = repmat(b, size(x));
    end
    if isscalar(c) 
        c = repmat(c, size(x));
    end
end

% Get indices for which the densty is defined
idx = find((x-c)>=0);

% Check size of the inputs
if (~isequal(size(x), size(a)) || ~isequal(size(x), size(b)) || ~isequal(size(x), size(c)))
    error('CheckInputArgs','All input arguments must have the same dimensions!')
end

if isempty(idx), return, end

% Compute the logged density

jdx = find( abs(a-1)<1e-12 & x>=c & (x-c)<1e-12) ;
ldens(jdx) = 1.0;
    
if ~isempty(idx)
    x0 = x(idx)-c(idx);
    x1 = x0./b(idx);
    x2 = x1.^a(idx);
    idx = setdiff(idx, jdx);
    ldens(idx) = log(a(idx)) - a(idx).*log(b(idx)) + (a(idx)-1).*log(x0) - x2 ;
end

% Compute the first and second derivatives.
if nargout>1
    x3 = (a(idx)-1)./x0;
    x4 = a(idx).*x2./x1./b(idx);
    Dldens(idx) = x3 - x4;
    if nargout>2
        D2ldens(idx) = -x3./x0 - (a(idx)-1).*x4./x1./b(idx);
    end
end

%@test:1
%$ try
%$    lpdfgweibull(1.0);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$    lpdfgweibull(1.0, .5);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ try
%$    lpdfgweibull(ones(2,2), .5*ones(2,1), .1*ones(3,1));
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ try
%$    a = lpdfgweibull(-1, .5, .1);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$ 
%$ if t(1)
%$    t(2) = isinf(a);
%$ end
%$
%$ T = all(t);
%@eof:4

%@test:5
%$ try
%$    a = lpdfgweibull(0, 1, 1);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(a-1.0)<1e-10;
%$ end
%$
%$ T = all(t);
%@eof:5

%@test:6
%$ try
%$    a = lpdfgweibull([0, -1], [1 1], [1 1]);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(a(1)-1.0)<1e-10;
%$    t(3) = isinf(a(2));
%$ end
%$
%$ T = all(t);
%@eof:6

%@test:7
%$ scale = 1;
%$ shape = 2;
%$ mode  = scale*((shape-1)/shape)^(1/shape);
%$ 
%$ try
%$    [a, b, c] = lpdfgweibull(mode, shape, scale);
%$    p = rand(1000,1)*4;
%$    am = lpdfgweibull(p, shape*ones(size(p)), scale*ones(size(p)));
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$
%$    t(2) = abs(b)<1e-8;
%$    t(3) = c<0;
%$    t(4) = all(am<a);
%$ end
%$
%$ T = all(t);
%@eof:7

%@test:8
%$ scale = 1;
%$ shape = 2;
%$ density  = @(x) exp(lpdfgweibull(x,shape,scale));
%$ 
%$ try
%$    if isoctave
%$        s = quadv(density, .0000000001, 100000,1e-10);
%$    else
%$        s = integral(density, 0, 100000);
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(s-1)<1e-6;
%$ end
%$
%$ T = all(t);
%@eof:8

%@test:9
%$ scale = 1;
%$ shape = 1;
%$ density  = @(x) exp(lpdfgweibull(x,shape,scale));
%$ 
%$ try
%$    if isoctave
%$        s = quadv(density, .0000000001, 100000,1e-10);
%$    else
%$        s = integral(density, 0, 100000);
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(s-1)<1e-6;
%$ end
%$
%$ T = all(t);
%@eof:9

%@test:10
%$ scale = 1;
%$ shape = .5;
%$ density  = @(x) exp(lpdfgweibull(x,shape,scale));
%$ 
%$ try
%$    if isoctave
%$        s = quadv(density, .0000000001, 100000,1e-10)
%$    else
%$        s = integral(density, 0, 100000);
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$        if isoctave()
%$            t(2) = abs(s-1)<5e-5;
%$        else
%$            t(2) = abs(s-1)<1e-6;
%$        end
%$ end
%$
%$ T = all(t);
%@eof:10

%@test:11
%$ scale = 1;
%$ shape = 2;
%$ xdens = @(x) x.*exp(lpdfgweibull(x,shape,scale));
%$ 
%$ try
%$    if isoctave
%$        s = quadv(xdens, .0000000001, 20,1e-10)
%$    else
%$        s = integral(xdens, 0, 100000);
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(s-scale*gamma(1+1/shape))<1e-6;
%$ end
%$
%$ T = all(t);
%@eof:11

%@test:12
%$ scale = 1;
%$ shape = 1;
%$ xdens = @(x) x.*exp(lpdfgweibull(x,shape,scale));
%$ 
%$ try
%$    if isoctave
%$        s = quadv(xdens, .0000000001, 100000,1e-10)
%$    else
%$        s = integral(xdens, 0, 100000);
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(s-scale*gamma(1+1/shape))<1e-6;
%$ end
%$
%$ T = all(t);
%@eof:12

%@test:13
%$ scale = 1;
%$ shape = .5;
%$ xdens = @(x) x.*exp(lpdfgweibull(x,shape,scale));
%$ 
%$ try
%$    if isoctave
%$        s = quadv(xdens, .0000000001, 100000,1e-10)
%$    else
%$        s = integral(xdens, 0, 100000);
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(s-scale*gamma(1+1/shape))<1e-6;
%$ end
%$
%$ T = all(t);
%@eof:13

%@test:14
%$ scale = 1;
%$ shape = 2;
%$ density = @(x) exp(lpdfgweibull(x,shape,scale));
%$ n = 200;
%$
%$ try
%$    s = NaN(n, 1);
%$    for i=1:n
%$        if isoctave()
%$            s(i) = quadv(density, .0000000001, .1*i,1e-10)
%$        else
%$            s(i) = integral(density, 0, .1*i);
%$        end
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    for i=1:n
%$        x = .1*i;
%$        q = 1-exp(-(x/scale)^shape);
%$        t(i+1) = abs(s(i)-q)<1e-6;
%$    end
%$ end
%$
%$ T = all(t);
%@eof:14

%@test:15
%$ scale = 1;
%$ shape = 1;
%$ density = @(x) exp(lpdfgweibull(x,shape,scale));
%$ n = 200;
%$
%$ try
%$    s = NaN(n, 1);
%$    for i=1:n
%$        if isoctave()
%$            s(i) = quadv(density, .0000000001, .1*i,1e-10)
%$        else
%$            s(i) = integral(density, 0, .1*i);
%$        end
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    for i=1:n
%$        x = .1*i;
%$        q = 1-exp(-(x/scale)^shape);
%$        t(i+1) = abs(s(i)-q)<1e-6;
%$    end
%$ end
%$
%$ T = all(t);
%@eof:15

%@test:16
%$ scale = 1;
%$ shape = .5;
%$ density = @(x) exp(lpdfgweibull(x,shape,scale));
%$ n = 200;
%$
%$ try
%$    s = NaN(n, 1);
%$    for i=1:n
%$        if isoctave()
%$            s(i) = quadv(density, .0000000001, .1*i,1e-10)
%$        else
%$            s(i) = integral(density, 0, .1*i);
%$        end
%$    end
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    for i=1:n
%$        x = .1*i;
%$        q = 1-exp(-(x/scale)^shape);
%$        if isoctave()
%$            t(i+1) = abs(s(i)-q)<5e-5;
%$        else
%$            t(i+1) = abs(s(i)-q)<1e-6;
%$        end
%$    end
%$ end
%$
%$ T = all(t);
%@eof:16