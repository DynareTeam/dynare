function B = append(A,a) % --*-- Unitary tests --*--

% append method for dates class.
%
% INPUTS 
%  o A    dates object.
%  o a    dates object with one element or string that can be interpreted as a date.
%
% OUTPUTS 
%  o B    dates object containing dates defined in A and a.
%
% EXAMPLE 1 
%  If A is a dates object with quarterly frequency, then B = A.append(dates('1950Q2')) and 
%  B = A.append('1950Q2') are equivalent syntaxes.

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

if isa(a,'dates')
    if ~isequal(length(a),1)
        error(['dates::append: Input argument ' inputname(2) ' has to be a dates object with one element.'])
    end
    if isempty(a)
        B = A;
        return
    end
elseif isdate(a)
    a = dates(a);
end

if ~isequal(A.freq, a.freq)
    error(['dates::append: A and a must have common frequency!'])
end

B = dates();
B.ndat = A.ndat+1;
B.freq = A.freq;
B.time = NaN(B.ndat,2);
B.time(1:A.ndat,:) = A.time;
B.time(A.ndat+1,:) = a.time; 

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009Q2';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine.
%$ d = dates(B4,B3,B2,B1);
%$ d = d.append(dates(B5));
%$
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009q2';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine.
%$ d = dates(B4,B3,B2,B1);
%$ d = d.append(B5);
%$
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:2