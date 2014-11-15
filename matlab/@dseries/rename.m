function ts = rename(ts,old,new) % --*-- Unitary tests --*--

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

if ~ischar(old) || ~ischar(new)
    error(['dseries::rename: Input arguments ''' inputname(2) ''' and ''' inputname(3) '''  have to be strings!'])
end
    
idname = find(strcmp(old,ts.name));

if isempty(idname)
    error(['dseries::rename: Variable ' old ' is unknown in dseries object ' inputname(1)  '!'])
end

ts.name(idname) = {new};

%@test:1
%$ t = zeros(8,1);
%$ ts = dseries([transpose(1:5), transpose(6:10)],'1950q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$ try
%$     ts = rename(ts,'Output','Production');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1950, 1]);
%$     t(5) = dassert(ts.vobs,2);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Production'; 'Consumption'});
%$     t(8) = dassert(ts.tex,{'Y_t'; 'C_t'});
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(8,1);
%$ ts = dseries([transpose(1:5), transpose(6:10)],'1950q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$ try
%$     ts = ts.rename('Output','Production');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1950, 1]);
%$     t(5) = dassert(ts.vobs,2);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Production'; 'Consumption'});
%$     t(8) = dassert(ts.tex,{'Y_t'; 'C_t'});
%$ end
%$
%$ T = all(t);
%@eof:2
