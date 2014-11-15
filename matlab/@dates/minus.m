function C = minus(A,B) % --*-- Unitary tests --*--

% Overloads the minus operator (-). If A and B are dates objects, the method returns the number of periods between A and B (so that A+C=B). If 
% one of the inputs is an integer or a vector of integers, the method shifts the dates object by X (the interger input) periods backward.

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
    if ~isequal(A.freq,B.freq)
        error(['dates::minus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must have common frequencies!'])
    end
    if isempty(A) || isempty(B)
        error(['dates::minus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must not be empty!'])
    end
    if ~isequal(length(A),length(B))
        if length(A)==1
            A.time = repmat(A.time,B.ndat,1);
            A.ndat = B.ndat;
        elseif length(B)==1
            B.time = repmat(B.time,A.ndat,1);
            B.ndat = A.ndat;
        else
            error(['dates::minus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' lengths are not consistent!'])
        end
    end
    C = zeros(length(A),1);
    id = find(~(A==B));
    if isempty(id)
        return
    end
    C(id) = A.time(id,2)-B.time(id,2) + (A.time(id,1)-B.time(id,1))*A.freq;
elseif (isvector(B) && isequal(length(B),A.ndat) && all(isint(B))) || isscalar(B) && isint(B) || isequal(length(A),1) && isvector(B) && all(isint(B))
    C = dates();
    C.freq = A.freq;
    C.time = add_periods_to_array_of_dates(A.time, A.freq, -B);
    C.ndat = rows(C.time);
else
    error('dates::plus: I don''t understand what you want to do! Check the manual.')
end

%@test:1
%$ % Define some dates objects
%$ d1 = dates('1950Q1','1950Q2','1960Q1');
%$ d2 = dates('1950Q3','1950Q4','1960Q1');
%$ d3 = dates('2000Q1');
%$ d4 = dates('2002Q2');
%$ % Call the tested routine.
%$ try
%$   e1 = d2-d1;
%$   e2 = d4-d3;
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,[2; 2; 0]);
%$   t(3) = dassert(e2,9);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates objects
%$ d1 = dates('1950Y','1951Y','1953Y');
%$ d2 = dates('1951Y','1952Y','1953Y');
%$ d3 = dates('2000Y');
%$ d4 = dates('1999Y');
%$ % Call the tested routine.
%$ try
%$   e1 = d2-d1;
%$   e2 = d4-d3;
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,[1; 1; 0]);
%$   t(3) = dassert(e2,-1);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dates objects
%$ d1 = dates('2000Y');
%$ d2 = dates('1999Y');
%$ % Call the tested routine.
%$ try
%$   e1 = d1-1;
%$   e2 = d2-(-1);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,d2);
%$   t(3) = dassert(e2,d1);
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define some dates objects
%$ d1 = dates('2000Q1');
%$ e1 = dates('1999Q4','1999Q3','1999Q2','1999Q1','1998Q4');
%$ % Call the tested routine.
%$ try
%$   f1 = d1-transpose(1:5);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,f1);
%$ end
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define some dates objects
%$ d1 = dates('1999Q4','1999Q3','1999Q2','1999Q1','1998Q4');
%$ e1 = dates('2000Q1')*5;
%$ % Call the tested routine.
%$ try
%$   f1 = d1-(-transpose(1:5));
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(e1,f1);
%$ end
%$ T = all(t);
%@eof:5