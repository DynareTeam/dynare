function A = set_names(B,varargin) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{A} =} times (@var{B},@code{varargin})
%! @anchor{@nSeries/set_names}
%! @sp 1
%! Specify names of the variables in a @ref{dseries} object.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item B
%! Dynare time series object instantiated by @ref{dseries}.
%! @item C
%! Dynare time series object instantiated by @ref{dseries}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object.
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

n = nargin-1;

if ~isdseries(B)
    error(['dseries::rename: ' inputname(1) ' must be a dseries object!'])
end

if ~isequal(vobs(B),n)
    error(['dseries::rename: The number of variables in ' inputname(1) ' does not match the number of declared names!'])
end

A = B;

for i=1:vobs(A)
    if ~isempty(varargin{i})
        A.name(i) = { varargin{i} };
    end
end

%@test:1
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Define names
%$ A_name = {'A1';'Variable_2';'A3'};
%$
%$ t = zeros(4,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],[],[]);
%$    ts2 = set_names(ts1,'A1',[],'A3');
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$ 
%$ if length(t)>1
%$    t(2) = dassert(ts2.vobs,3);
%$    t(3) = dassert(ts2.nobs,10);
%$    t(4) = dassert(ts2.name,A_name);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Define names
%$ A_name = {'A1';'Variable_2';'A3'};
%$
%$ t = zeros(4,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],[],[]);
%$    ts1 = ts1.set_names('A1',[],'A3');
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name,A_name);
%$ end
%$ T = all(t);
%@eof:2
