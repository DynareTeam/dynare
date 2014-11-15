function C = plus(A,B) % --*-- Unitary tests --*--

% Overloads the plus operator. If A and B are dates objects the method combines A and B without removing repetitions. If
% one of the inputs is an integer or a vector of integers, the method shifts the dates object by X (the interger input)
% periods forward.

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

if isa(B,'dates')
    % Concatenate dates objects without removing repetitions if A and B are not disjoint sets of dates.
    if ~isequal(A.freq,B.freq) && A.ndat>0 && B.ndat>0
        error(['dates::plus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must have common frequencies!'])
    end
    if isempty(B)
        C = A;
        return
    end
    if isempty(A)
        C = B;
        return
    end
    C = dates();
    C.freq = A.freq;
    C.time = [A.time; B.time];
    C.ndat = A.ndat+B.ndat;
elseif (isvector(B) && isequal(length(B),A.ndat) && all(isint(B))) || isscalar(B) && isint(B) || isequal(length(A),1) && isvector(B) && all(isint(B))
    C = dates();
    C.freq = A.freq;
    C.time = add_periods_to_array_of_dates(A.time, A.freq, B);
    C.ndat = rows(C.time);
else
    error('dates::plus: I don''t understand what you want to do! Check the manual.')
end

%@test:1
%$ % Define some dates objects
%$ d1 = dates('1950Q1','1950Q2') ;
%$ d2 = dates('1950Q3','1950Q4') ;
%$ d3 = dates('1950Q1','1950Q2','1950Q3','1950Q4') ;
%$
%$ % Call the tested routine.
%$ try
%$   e1 = d1+d2;
%$   e2 = d1+d2+d3;
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1) 
%$   t(2) = dassert(e1,d3);
%$   t(3) = dassert(e2,dates('1950Q1','1950Q2','1950Q3','1950Q4','1950Q1','1950Q2','1950Q3','1950Q4'));
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates objects
%$ d1 = dates('1950Q1');
%$ e1 = dates('1950Q2');
%$ e2 = dates('1950Q3');
%$ e3 = dates('1950Q4');
%$ e4 = dates('1951Q1');
%$ e5 = dates('1950Q2','1950Q3','1950Q4','1951Q1');
%$
%$ % Call the tested routine.
%$ try
%$   f1 = d1+1;
%$   f2 = d1+2;
%$   f3 = d1+3;
%$   f4 = d1+4;
%$   f5 = d1+transpose(1:4);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,f1);
%$   t(3) = dassert(e2,f2);
%$   t(4) = dassert(e3,f3);
%$   t(5) = dassert(e4,f4);
%$   t(6) = dassert(e5,f5);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dates objects
%$ d1 = dates('1950Q1');
%$ e1 = dates('1949Q4');
%$ e2 = dates('1949Q3');
%$ e3 = dates('1949Q2');
%$ e4 = dates('1949Q1');
%$ e5 = dates('1948Q4');
%$
%$ % Call the tested routine.
%$ try
%$   f1 = d1+(-1);
%$   f2 = d1+(-2);
%$   f3 = d1+(-3);
%$   f4 = d1+(-4);
%$   f5 = d1+(-5);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,f1);
%$   t(3) = dassert(e2,f2);
%$   t(4) = dassert(e3,f3);
%$   t(5) = dassert(e4,f4);
%$   t(6) = dassert(e5,f5);
%$ end
%$ T = all(t);
%@eof:3