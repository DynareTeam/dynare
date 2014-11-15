function C = ge(A,B)  % --*-- Unitary tests --*--

% Overloads the >= operator for dates objects.
%
% INPUTS 
%  o A    dates object with n or 1 elements.
%  o B    dates object with n or 1 elements.
%
% OUTPUTS 
%  o C    column vector of max(n,1) elements (zeros or ones).

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

if ~isequal(nargin,2)
    error('dates::ge: I need exactly two input arguments!')
end

if ~isa(A,'dates') || ~isa(B,'dates')
    error(['dates::ge: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' have to be a dates objects!'])
end

if ~isequal(A.freq,B.freq)
    C = false;
    return
end

if isequal(A.ndat, B.ndat)
    C = (A==B);
    idx = find(C==0);
    for i=1:length(idx)
        C(idx(i)) = greaterorequal(A.time(idx(i),:), B.time(idx(i),:));
    end
else
    if isequal(A.ndat,1)
        C = false(B.ndat,1);
        for i=1:B.ndat
            C(i) = greaterorequal(A.time, B.time(i,:));
        end
    elseif isequal(B.ndat,1)
        C = false(A.ndat,1);
        for i=1:A.ndat
            C(i) = greaterorequal(A.time(i,:), B.time);
        end
    else
        C = false;
    end
end


function c = greaterorequal(a,b)
    if a(1)>b(1)
        c = true;
    else
        if a(1)<b(1)
            c = false;
        else
            if a(2)>=b(2)
                c = true;
            else
                c = false;
            end
        end
    end


%@test:1
%$ % Define some dates
%$ date_2 = '1950Q2';
%$ date_3 = '1950Q3';
%$ date_4 = '1950Q1';
%$ date_5 = '1949Q2';
%$
%$ % Call the tested routine.
%$ d2 = dates(date_2);
%$ d3 = dates(date_3);
%$ d4 = dates(date_4);
%$ d5 = dates(date_5);
%$ i1 = (d2>=d3);
%$ i2 = (d3>=d4);
%$ i3 = (d4>=d2);
%$ i4 = (d5>=d4);
%$ i5 = (d5>=d5); 
%$
%$ % Check the results.
%$ t(1) = dassert(i1,false);
%$ t(2) = dassert(i2,true);
%$ t(3) = dassert(i3,false);
%$ t(4) = dassert(i4,false);
%$ t(5) = dassert(i5,true);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ B1 = '1945Q1';
%$ B2 = '1945Q2';
%$ B3 = '1945Q3';
%$ B4 = '1945Q4';
%$ B5 = '1950Q1';
%$
%$ % Create dates objects.
%$ dd = dates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dassert(dates(B1)>=dates(B2),false);
%$ t(2) = dassert(dates(B2)>=dates(B1),true);
%$ t(3) = dassert(dates(B2)>=dates(B2),true);
%$ t(4) = dassert(dd>=dates(B5),false(4,1));
%$ t(5) = dassert(dates(B5)>=dd,true(4,1));
%$ t(6) = dassert(dates(B1)>=dd,[true; false(3,1)]);
%$ T = all(t);
%@eof:2