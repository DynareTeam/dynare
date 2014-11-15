function us = ydiff(ts) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{us} =} ydiff (@var{ts})
%! @anchor{ydiff}
%! @sp 1
%! Computes annual differences.
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

us = ts;

switch frequency(ts)
  case 1
    us.data(2:end,:) = ts.data(2:end,:)-ts.data(1:end-1,:);
    us.data(1,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ydiff(' us.name{i} ')']};
        us.tex(i) = {['\Delta ' us.tex{i}]};
    end
  case 4
    us.data(5:end,:) = ts.data(5:end,:)-ts.data(1:end-4,:);
    us.data(1:4,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ydiff(' us.name{i} ')']};
        us.tex(i) = {['\Delta_4 ' us.tex{i}]};
    end
  case 12
    us.data(13:end,:) = ts.data(13:end,:)-ts.data(1:end-12,:);
    us.data(1:12,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ydiff(' us.name{i} ')']};
        us.tex(i) = {['\Delta_{12} ' us.tex{i}]};
    end
  case 52
    us.data(53:end,:) = ts.data(53:end,:)-ts.data(1:end-52,:);
    us.data(1:52,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ydiff(' us.name{i} ')']};
        us.tex(i) = {['\Delta_{52} ' us.tex{i}]};
    end
  otherwise
    error(['dseries::ygrowth: object ' inputname(1) ' has unknown frequency']);
end

%@test:1
%$ t = zeros(4,1);
%$
%$ try
%$     data = transpose(1:100);
%$     ts = dseries(data,'1950Q1',{'A1'},{'A_1'});
%$     ts = ts.ydiff;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(4,ts.vobs);
%$     DATA = [DATA; 4*ones(ts.nobs-4,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$     t(3) = dassert(ts.name{1},['ydiff(A1)']);
%$     t(4) = dassert(ts.tex{1},['\\Delta_4 A_1']);
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(4,1);
%$
%$ try
%$     data = transpose(1:100);
%$     ts = dseries(data,'1950M1',{'A1'},{'A_1'});
%$     ts = ts.ydiff;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(12,ts.vobs);
%$     DATA = [DATA; 12*ones(ts.nobs-12,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$     t(3) = dassert(ts.name{1},['ydiff(A1)']);
%$     t(4) = dassert(ts.tex{1},['\\Delta_{12} A_1']);
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(2,1);
%$
%$ try
%$     data = transpose(1:100);
%$     ts = dseries(data,'1950W1',{'A1'},{'A_1'});
%$     ts = ts.ydiff;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(52,ts.vobs);
%$     DATA = [DATA; 52*ones(ts.nobs-52,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$     t(3) = dassert(ts.name{1},['ydiff(A1)']);
%$     t(4) = dassert(ts.tex{1},['\\Delta_{52} A_1']);
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ t = zeros(4,1);
%$
%$ try
%$     data = transpose(1:100);
%$     ts = dseries(data,'1950Y',{'A1'},{'A_1'});
%$     ts = ts.ydiff;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(1,ts.vobs);
%$     DATA = [DATA; ones(ts.nobs-1,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$     t(3) = dassert(ts.name{1},['ydiff(A1)']);
%$     t(4) = dassert(ts.tex{1},['\\Delta A_1']);
%$ end
%$
%$ T = all(t);
%@eof:4
