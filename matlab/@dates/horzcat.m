function B = horzcat(varargin) % --*-- Unitary tests --*--

% Overloads the horzcat method for dates objects.
%
% INPUTS 
%  o A1    dates object.
%  o A2    dates object.
%  o ...
%
% OUTPUTS 
%  o B    dates object containing dates defined in A1, A2, ...
%
% EXAMPLE 1 
%  If A, B and C are dates objects the following syntax:
%    
%    D = [A, B, C] ;
%
%  Defines a dates object D containing the dates appearing in A, B and C.
    
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

if ~all(cellfun(@isdates,varargin))
    error('dates::horzcat: All input arguments must be dates objects.')
end

n = nargin;
B = varargin{1};

if isequal(n,1), return, end

for i=2:n
    C = varargin{i};
    if isequal(B.freq,C.freq)
        if ~isempty(C)
            B.ndat = B.ndat + C.ndat;
            B.time = [B.time; C.time];
        end
    else
        error('dates::horzcat: All input arguments must have the same frequency!')
    end
end

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
%$ d = [d, dates(B5)];
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
%$ B5 = '2009Q2';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine.
%$ d = dates(B4,B3,B2);
%$ d = [d, dates(B1), dates(B5)];
%$
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:2

%@test:3
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
%$ d = dates(B4,B3,B2);
%$ d = [d, dates(B1,B5)];
%$
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:3

%@test:4
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
%$ d = dates(B4,B3,B2);
%$ d = [d, [dates(B1), dates(B5)]];
%$
%$ % Check the results.
%$ t(1) = dassert(d.time,e.time);
%$ t(2) = dassert(d.freq,e.freq);
%$ t(3) = dassert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:4