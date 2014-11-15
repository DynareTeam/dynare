function [a,b] = align(a, b) % --*-- Unitary tests --*--
    
%@info:
%! @deftypefn {Function File} {[@var{a}, @var{b}] =} align (@var{a}, @var{b})
%! @anchor{dseries/align}
%! @sp 1
%! If dseries objects @var{a} and @var{b} are defined on different time ranges, extend @var{a} and/or
%! @var{b} with NaNs so that they are defined on the same time range.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Object instantiated by @ref{dseries}.
%! @item b
%! Object instantiated by @ref{dseries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Object instantiated by @ref{dseries}.
%! @item b
%! Object instantiated by @ref{dseries}.
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2013 Dynare Team
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

if ~isequal(frequency(a),frequency(b))
    error(['dseries::align: ''' inputname(1) ''' and ''' inputname(2) ''' dseries objects must have common frequencies!'])
end

init = min(firstdate(a),firstdate(b));
last = max(lastdate(a),lastdate(b));

if isempty(intersect(a.dates,b.dates))
    error(['dseries::align: ''' inputname(1) ''' and ''' inputname(2) ''' dseries object must have at least one common date!'])
end

a_init = init;
b_init = init;
a_last = last;
b_last = last;

if firstdate(b)>init
    n = firstdate(b)-init;
    b.data = [NaN(n, vobs(b)); b.data];
    b_init = init;
end

if firstdate(a)>init
    n = firstdate(a)-init;
    a.data = [NaN(n, vobs(a)); a.data];
    a_init = init;
end

if lastdate(b)<last
    n = last-lastdate(b);
    b.data = [b.data; NaN(n, vobs(b))];
end

if lastdate(a)<last
    n = last-lastdate(a);
    a.data = [a.data; NaN(n, vobs(a))];
end

a.dates = a_init:a_init+(nobs(a)-1);
b.dates = b_init:b_init+(nobs(b)-1);

%@test:1
%$ % Define a datasets.
%$ A = rand(8,3); B = rand(7,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial dates
%$ A_init = '1990Q1';
%$ B_init = '1989Q2';
%$
%$ % Instantiate two dseries objects
%$ ts1 = dseries(A,A_init,A_name);
%$ ts2 = dseries(B,B_init,B_name);
%$
%$ try
%$   [ts1, ts2] = align(ts1, ts2);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(ts1.nobs,ts2.nobs);
%$   t(3) = dassert(ts1.init,ts2.init);
%$   t(4) = dassert(ts1.data,[NaN(3,3); A], 1e-15);
%$   t(5) = dassert(ts2.data,[B; NaN(4,2)], 1e-15);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(8,3); B = rand(7,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial dates
%$ A_init = '1990Q1';
%$ B_init = '1990Q1';
%$
%$ % Instantiate two dseries objects
%$ ts1 = dseries(A,A_init,A_name);
%$ ts2 = dseries(B,B_init,B_name);
%$
%$ try
%$   [ts1, ts2] = align(ts1, ts2);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(ts1.nobs,ts2.nobs);
%$   t(3) = dassert(ts1.init,ts2.init);
%$   t(4) = dassert(ts1.data,A, 1e-15);
%$   t(5) = dassert(ts2.data,[B; NaN(1,2)], 1e-15);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a datasets.
%$ A = rand(8,3); B = rand(7,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial dates
%$ A_init = '1990Q1';
%$ B_init = '1990Q1';
%$
%$ % Instantiate two dseries objects
%$ ts1 = dseries(A,A_init,A_name);
%$ ts2 = dseries(B,B_init,B_name);
%$
%$ try
%$   [ts2, ts1] = align(ts2, ts1);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(ts1.nobs,ts2.nobs);
%$   t(3) = dassert(ts1.init,ts2.init);
%$   t(4) = dassert(ts1.data,A, 1e-15);
%$   t(5) = dassert(ts2.data,[B; NaN(1,2)], 1e-15);
%$ end
%$ T = all(t);
%@eof:3