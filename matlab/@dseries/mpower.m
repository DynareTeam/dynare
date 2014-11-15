function A = mpower(B,C) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{A} =} mpower (@var{B},@var{C})
%! @anchor{@dseries/mpower}
%! @sp 1
%! Overloads the mpower method for the Dynare time series class (@ref{dseries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item B
%! Dynare time series object instantiated by @ref{dseries}, with T observations and N variables.
%! @item C
%! Real scalar or a dseries object with T observations and N variables.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item A
%! dseries object with T observations and N variables.
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

if isnumeric(B) && isvector(B) && length(B)>1
    if ~isdseries(C)
        error('dseries::mpower: Second input argument must be a dseries object!')
    end
    A = C;
    A.data = bsxfun(@power,C.data,B);
    return;
end

if isnumeric(C) && isvector(C) && length(C)>1
    if ~isdseries(B)
        error('dseries::mpower: First input argument must be a dseries object!')
    end
    A = B;
    A.data = bsxfun(@power,B.data,C);
    return
end

if isdseries(B) && isnumeric(C) && isreal(C) &&  isscalar(C)
    A = dseries();
    A.dates = B.dates;
    A.data = B.data.^C;
    A.name = cell(vobs(A),1);
    A.tex = cell(vobs(A),1);
    for i=1:vobs(A)
        A.name(i) = {['power(' B.name{i} ';' num2str(C) ')']};
        A.tex(i) = {[B.tex{i} '^' num2str(C) ]};
    end
    return
end

if isdseries(B) && isdseries(C)
    if isequal(nobs(B),nobs(C)) && isequal(vobs(B), vobs(C)) && isequal(frequency(B),frequency(C))
        A = dseries();
        A.data = B.data.^C.data;
        A.dates = B.dates;
        A.name = cell(vobs(A),1);
        A.tex = cell(vobs(A),1);
        for i=1:vobs(A)
            A.name(i) = {['power(' B.name{i} ';' C.name{i} ')']};
            A.tex(i) = {[B.tex{i} '^{' C.tex{i} '}']};
        end
    else
        error('dseries::mpower: If both input arguments are dseries objects, they must have the same numbers of variables and observations and common frequency!')
    end
    return
end

error(['dseries::mpower: Wrong calling sequence!'])

%@test:1
%$ % Define a datasets.
%$ A = rand(10,2); B = randn(10,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2'}; B_name = {'B1';'B2'};
%$
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = dseries(B,[],B_name,[]);
%$    ts3 = ts1^ts2;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts3.vobs,2);
%$    t(3) = dassert(ts3.nobs,10);
%$    t(4) = dassert(ts3.data,A.^B,1e-15);
%$    t(5) = dassert(ts3.name,{'power(A1;B1)';'power(A2;B2)'});
%$    t(6) = dassert(ts3.tex,{'A1^{B1}';'A2^{B2}'});
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts3 = ts1^2;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts3.vobs,2);
%$    t(3) = dassert(ts3.nobs,10);
%$    t(4) = dassert(ts3.data,A.^2,1e-15);
%$    t(5) = dassert(ts3.name,{'power(A1;2)';'power(A2;2)'});
%$    t(6) = dassert(ts3.tex,{'A1^2';'A2^2'});
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a dseries object
%$ ts1=dseries([1 1;2 2;3 3], '1999y', {'MyVar1','MyVar2'});
%$
%$ % Use the power
%$ try
%$    ts2 = ts1^transpose(1:3);
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts2.vobs,2);
%$    t(3) = dassert(ts2.nobs,3);
%$    t(4) = dassert(ts2.data,bsxfun(@power,ts1.data,transpose(1:3)),1e-15);
%$    t(5) = dassert(ts2.name,{'MyVar1';'MyVar2'});
%$    t(6) = dassert(ts2.tex,{'MyVar1';'MyVar2'});
%$ end
%$ T = all(t);
%@eof:3