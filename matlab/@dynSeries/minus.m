function A = minus(B,C)
%@info:
%! @deftypefn {Function File} {@var{A} =} minus (@var{B},@var{C})
%! @anchor{@dynSeries/minus}
%! @sp 1
%! Overloads the minus method for the Dynare time series class (@ref{dynSeries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item B
%! Dynare time series object instantiated by @ref{dynSeries}.
%! @item C
%! Dynare time series object instantiated by @ref{dynSeries}.
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
[B, C] = align(B, C);
if ~isequal(B.vobs,C.vobs) && ~(isequal(B.vobs,1) || isequal(C.vobs,1))
    error(['dynSeries::minus: Cannot substract ' inputname(1) ' and ' inputname(2) ' (wrong number of variables)!'])
else
    if B.vobs>C.vobs
        idB = 1:B.vobs;
        idC = ones(1:B.vobs);
    elseif B.vobs<C.vobs
        idB = ones(1,C.vobs);
        idC = 1:C.vobs;
    else
        idB = 1:B.vobs;
        idC = 1:C.vobs;
    end
end

if ~isequal(B.freq,C.freq)
    error(['dynSeries::plus: Cannot substract ' inputname(1) ' and ' inputname(2) ' (frequencies are different)!'])
end

if ~isequal(B.nobs,C.nobs) || ~isequal(B.init,C.init)
    [B, C] = align(B, C);
end

if isempty(B)
    A = -C;
    return
end

if isempty(C)
    A = B;
    return
end

A = dynSeries();

A.freq = B.freq;
A.init = B.init;
A.time = B.time;
A.nobs = max(B.nobs,C.nobs);
A.vobs = max(B.vobs,C.vobs);
A.name = cell(A.vobs,1);
A.tex = cell(A.vobs,1);
for i=1:A.vobs
    A.name(i) = {['minus(' B.name{idB(i)} ',' C.name{idC(i)} ')']};
    A.tex(i) = {['(' B.tex{idB(i)} '-' C.tex{idC(i)} ')']};
end
A.data = bsxfun(@minus,B.data,C.data);

%@test:1
%$ % Define a datasets.
%$ A = rand(10,2); B = randn(10,1);
%$
%$ % Define names
%$ A_name = {'A1';'A2'}; B_name = {'B1'};
%$
%$ t = zeros(5,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    ts2 = dynSeries(B,[],B_name,[]);
%$    ts3 = ts1-ts2;
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(ts3.vobs,2);
%$    t(3) = dyn_assert(ts3.nobs,10);
%$    t(4) = dyn_assert(ts3.data,[A(:,1)-B, A(:,2)-B],1e-15);
%$    t(5) = dyn_assert(ts3.name,{'minus(A1,B1)';'minus(A2,B1)'});
%$ end
%$ T = all(t);
%@eof:1

%@test:3
%$ % Define a datasets.
%$ A = rand(10,2); B = randn(5,1);
%$
%$ % Define names
%$ A_name = {'A1';'A2'}; B_name = {'B1'};
%$
%$ t = zeros(5,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    ts2 = dynSeries(B,[],B_name,[]);
%$    ts3 = ts1-ts2;
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(ts3.vobs,2);
%$    t(3) = dyn_assert(ts3.nobs,10);
%$    t(4) = dyn_assert(ts3.data,[A(1:5,1)-B(1:5), A(1:5,2)-B(1:5) ; NaN(5,2)],1e-15);
%$    t(5) = dyn_assert(ts3.name,{'minus(A1,B1)';'minus(A2,B1)'});
%$ end
%$ T = all(t);
%@eof:3
