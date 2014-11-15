function us = ygrowth(ts) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{us} =} ygrowth (@var{ts})
%! @anchor{ygrowth}
%! @sp 1
%! Computes annual growth rates.
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
    us.data(2:end,:) = ts.data(2:end,:)./ts.data(1:end-1,:) - 1;
    us.data(1,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ygrowth(' us.name{i} ')']};
        us.tex(i) = {['\delta ' us.tex{i}]};
    end
  case 4
    us.data(5:end,:) = ts.data(5:end,:)./ts.data(1:end-4,:) - 1;
    us.data(1:4,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ygrowth(' us.name{i} ')']};
        us.tex(i) = {['\delta_4 ' us.tex{i}]};
    end
  case 12
    us.data(13:end,:) = ts.data(13:end,:)./ts.data(1:end-12,:) - 1;
    us.data(1:12,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ygrowth(' us.name{i} ')']};
        us.tex(i) = {['\delta_{12} ' us.tex{i}]};
    end
  case 52
    us.data(53:end,:) = ts.data(53:end,:)./ts.data(1:end-52,:) - 1;
    us.data(1:52,:) = NaN;
    for i = 1:vobs(ts)
        us.name(i) = {['ygrowth(' us.name{i} ')']};
        us.tex(i) = {['\delta_{52} ' us.tex{i}]};
    end
  otherwise
    error(['dseries::ygrowth: object ' inputname(1) ' has unknown frequency']);
end

%@test:1
%$ t = zeros(2,1);
%$
%$ try
%$     data = repmat(transpose(1:4),100,1);
%$     ts = dseries(data,'1950Q1');
%$     ts = ts.ygrowth;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(4,ts.vobs);
%$     DATA = [DATA; zeros(ts.nobs-4,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(2,1);
%$
%$ try
%$     data = repmat(transpose(1:12),100,1);
%$     ts = dseries(data,'1950M1');
%$     ts = ts.ygrowth;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(12,ts.vobs);
%$     DATA = [DATA; zeros(ts.nobs-12,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(2,1);
%$
%$ try
%$     data = repmat(transpose(1:52),100,1);
%$     ts = dseries(data,'1950W1');
%$     ts = ts.ygrowth;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$
%$ if length(t)>1
%$     DATA = NaN(52,ts.vobs);
%$     DATA = [DATA; zeros(ts.nobs-52,ts.vobs)];
%$     t(2) = dassert(ts.data,DATA);
%$ end
%$
%$ T = all(t);
%@eof:3
