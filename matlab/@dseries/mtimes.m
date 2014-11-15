function A = mtimes(B,C) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{A} =} mtimes (@var{B},@var{C})
%! @anchor{@dseries/mtimes}
%! @sp 1
%! Overloads the mtimes method for the Dynare time series class (@ref{dseries}).
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

% Copyright (C) 2012-2014 Dynare Team
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

if isnumeric(B) && (isscalar(B) ||  isvector(B))
    if ~isdseries(C)
        error('dseries::mtimes: Second input argument must be a dseries object!')
    end
    A = C;
    A.data = bsxfun(@times,C.data,B);
    return;
end

if isnumeric(C) && (isscalar(C) || isvector(C))
    if ~isdseries(B)
        error('dseries::mtimes: First input argument must be a dseries object!')
    end
    A = B;
    A.data = bsxfun(@times,B.data,C);
    return
end

if isdseries(B) && isdseries(C)
    % Element by element multiplication of two dseries object
    if ~isequal(vobs(B), vobs(C)) && ~(isequal(vobs(B),1) || isequal(vobs(C),1))
        error(['dseries::times: Cannot multiply ' inputname(1) ' and ' inputname(2) ' (wrong number of variables)!'])
    else
        if vobs(B)>vobs(C)
            idB = 1:vobs(B);
            idC = ones(1:vobs(B));
        elseif vobs(B)<vobs(C)
            idB = ones(1,vobs(C));
            idC = 1:vobs(C);
        else
            idB = 1:vobs(B);
            idC = 1:vobs(C);
        end
    end
    if ~isequal(frequency(B),frequency(C))
        error(['dseries::times: Cannot multiply ' inputname(1) ' and ' inputname(2) ' (frequencies are different)!'])
    end
    if ~isequal(nobs(B), nobs(C)) || ~isequal(firstdate(B),firstdate(C))
        [B, C] = align(B, C);
    end
    A = dseries();
    A.dates = B.dates;
    A_vobs = max(vobs(B),vobs(C));
    A.name = cell(A_vobs,1);
    A.tex = cell(A_vobs,1);
    for i=1:A_vobs
        A.name(i) = {['multiply(' B.name{idB(i)} ';' C.name{idC(i)} ')']};
        A.tex(i) = {['(' B.tex{idB(i)} '*' C.tex{idC(i)} ')']};
    end
    A.data = bsxfun(@times,B.data,C.data);
else
    error()
end

%@test:1
%$ % Define a datasets.
%$ A = rand(10,2); B = randn(10,1);
%$
%$ % Define names
%$ A_name = {'A1';'A2'}; B_name = {'B1'};
%$
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = dseries(B,[],B_name,[]);
%$    ts3 = ts1*ts2;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts3.vobs,2);
%$    t(3) = dassert(ts3.nobs,10);
%$    t(4) = dassert(ts3.data,[A(:,1).*B, A(:,2).*B],1e-15);
%$    t(5) = dassert(ts3.name,{'multiply(A1;B1)';'multiply(A2;B1)'});
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,2); B = pi;
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = ts1*B;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts2.vobs,2);
%$    t(3) = dassert(ts2.nobs,10);
%$    t(4) = dassert(ts2.data,A*B,1e-15);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a datasets.
%$ A = rand(10,2); B = pi;
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = B*ts1;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts2.vobs,2);
%$    t(3) = dassert(ts2.nobs,10);
%$    t(4) = dassert(ts2.data,A*B,1e-15);
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a datasets.
%$ A = rand(10,2); B = A(1,:);
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = B*ts1;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = dassert(ts2.vobs,2);
%$    t(3) = dassert(ts2.nobs,10);
%$    t(4) = dassert(ts2.data,bsxfun(@times,A,B),1e-15);
%$ end
%$ T = all(t);
%@eof:4
