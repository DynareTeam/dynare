function vs = chain(ts,us)  % --*-- Unitary tests --*--

% Copyright (C) 2014 Dynare Team
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

if vobs(ts)-vobs(us)
    error(['dseries::chain: dseries objects ' inputname(1) ' and ' inputname(2) ' must have the same number of variables!'])
end

if frequency(ts)-frequency(us)
    error(['dseries::chain: dseries objects ' inputname(1) ' and ' inputname(2) ' must have common frequencies!'])
end

if lastdate(ts)<firstdate(us)
    error(['dseries::chain: The last date in ' inputname(1) ' (' date2string(ts.dates(end)) ') must not preceed the first date in ' inputname(2) ' (' date2string(us.dates(1)) ')!'])
end

tdx = find(sum(bsxfun(@eq,us.dates.time,ts.dates.time(end,:)),2)==2);
GrowthFactor = us.data(tdx+1:end,:)./us.data(tdx:end-1,:);
CumulatedGrowthFactors = cumprod(GrowthFactor);

vs = ts;
vs.data = [vs.data; bsxfun(@times,CumulatedGrowthFactors,vs.data(end,:))];

vs.dates = firstdate(vs):firstdate(vs)+nobs(vs);

%@test:1
%$ try
%$     ts = dseries([1; 2; 3; 4],dates('1950Q1')) ;
%$     us = dseries([3; 4; 5; 6],dates('1950Q3')) ;
%$     vs = chain(ts,us);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(vs.freq,4);
%$     t(3) = dassert(vs.init.freq,4);
%$     t(4) = dassert(vs.init.time,[1950, 1]);
%$     t(5) = dassert(ts.vobs,1);
%$     t(6) = dassert(vs.nobs,6);
%$     t(7) = isequal(vs.data,transpose(1:6));
%$ end
%$
%$ T = all(t);
%@eof:1