function us = lag(ts,p) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{us} =} lag (@var{ts})
%! @anchor{lag}
%! @sp 1
%! Computes lagged time series.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @var
%! @item ts
%! Dynare time series object, instantiated by @ref{dseries}
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item us
%! Dynare time series object with transformed data field.
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

% Set default number of lags
if nargin<2
    p = 1;
end

if p<=0
    error('dseries::lag: Second input argument must be strictly positive! Use lead method instead.')
end

% Copy of ts dseries object
us = ts;

% Update data member
us.data = [NaN(p, vobs(ts));  ts.data(1:end-p,:)];

for i=1:vobs(ts)
    us.name(i) = {[ 'lag(' ts.name{i} ',' int2str(p) ')']};
    us.tex(i) = {[ ts.tex{i} '_{-' int2str(p) '}']};
end

%@test:1
%$ t = zeros(4,1);
%$
%$ try
%$     data = transpose(0:1:50);
%$     ts = dseries(data,'1950Q1');
%$     a = ts.lag;
%$     b = ts.lag.lag;
%$     c = lag(ts,2);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     DATA = [NaN(1,ts.vobs); transpose(0:1:49)];
%$     t(2) = dassert(a.data,DATA,1e-15);
%$     DATA = [NaN(2,ts.vobs); transpose(0:1:48)];
%$     t(3) = dassert(b.data,DATA,1e-15);
%$     t(4) = dassert(b.data,c.data,1e-15);
%$ end
%$
%$ T = all(t);
%@eof:1