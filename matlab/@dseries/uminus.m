function A = uminus(B) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{A} =} plus (@var{B},@var{C})
%! @anchor{@dseries/uminus}
%! @sp 1
%! Overloads the uminus method for the Dynare time series class (@ref{dseries}).
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

% Copyright (C) 2012-2013 Dynare Team
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

A = dseries();

A.data = -(B.data);
A.dates = B.dates;

A.name = cell(vobs(A),1);
A.tex = cell(vobs(A),1);
for i = 1:vobs(A)
    A.name(i) = {[ '-' B.name{i}]};
    A.tex(i) = {[ '-' B.tex{i}]};
end

%@test:1
%$ % Define a datasets.
%$ A = rand(10,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ A_tex = {'A_1';'A_2'};
%$ t = zeros(6,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,A_tex);
%$    ts2 = -ts1;
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(ts2.vobs,2);
%$    t(3) = dassert(ts2.nobs,10);
%$    t(4) = dassert(ts2.data,-A,1e-15);
%$    t(5) = dassert(ts2.name,{'-A1';'-A2'});
%$    t(6) = dassert(ts2.tex,{'-A_1';'-A_2'});
%$ end
%$ T = all(t);
%@eof:1
