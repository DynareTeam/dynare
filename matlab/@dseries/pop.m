function [ts,id] = pop(ts,a) % --*-- Unitary tests --*--

% Removes a variable from a dseries object.

%@info:
%! @deftypefn {Function File} {@var{ts} =} pop (@var{ts}, @var{a})
%! @anchor{dseries/pop}
%! @sp 1
%! Pop method for the dseries class. Removes a variable from a dseries object.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Object instantiated by @ref{dseries}.
%! @item a
%! String, name of the variable to be removed.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Object instantiated by @ref{dseries}, without variable (@var{a}).
%! @item id
%! Scalar integer, position of variable (@var{a}) in the original dseries object @var{ts}.
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

if nargin<2
    % Removes the last variable
    id = vobs(ts);
else
    id = find(strcmp(a,ts.name));
end

if isempty(id)
    id = 0;
    return
end

ts.data(:,id) = [];
ts.name(id) = [];
ts.tex(id) = [];

%@test:1
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$
%$ t = zeros(4,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = pop(ts1,'A2');
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(ts2.vobs,2);
%$    t(3) = dassert(ts2.nobs,10);
%$    t(4) = dassert(ts2.data,[A(:,1), A(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$
%$ t = zeros(2,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    [ts2,id] = pop(ts1,'A4');
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(id,0);
%$    t(2) = dassert(ts1,ts2);
%$ end
%$ T = all(t);
%@eof:2
