function dd = sort(dd) % --*-- Unitary tests --*--

% Sort method for dates class.
%
% INPUTS 
%  o dd    dates object.
%
% OUTPUTS 
%  o dd    dates object (with dates sorted by increasing order).

% Copyright (C) 2011-2013 Dynare Team
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

if isequal(dd.ndat,1)
    return
end

dd.time = sortrows(dd.time,[1,2]);

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4];
%$ e.freq = 4;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dates(B1,B2,B3,B4);
%$ d = d.sort;
%$ 
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4];
%$ e.freq = 4;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dates(B1,B2,B3,B4);
%$ d = d.sort();
%$
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1