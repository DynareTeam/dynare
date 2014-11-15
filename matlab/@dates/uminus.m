function B = uminus(A)  % --*-- Unitary tests --*--

% Overloads the unary minus operator for dates objects. Shifts all the elements by one period.
%
% INPUTS 
%  o A    dates object with n elements.
%
% OUTPUTS 
%  o B    dates object with n elements.

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

B = dates(A);
B.time(:,2) = B.time(:,2)-1;
idx = find(B.time(:,2)==0);
B.time(idx,1) = B.time(idx,1)-1;
B.time(idx,2) = B.freq;

%@test:1
%$ % Define some dates
%$ date_1 = '1950Y';
%$ date_2 = '1950Q2';
%$ date_3 = '1950Q1';
%$ date_4 = '1950M2';
%$ date_5 = '1950M1';
%$
%$ % Call the tested routine.
%$ d1 = dates(date_1); d1 = -d1;
%$ d2 = dates(date_2); d2 = -d2;
%$ d3 = dates(date_3); d3 = -d3;
%$ d4 = dates(date_4); d4 = -d4;
%$ d5 = dates(date_5); d5 = -d5;
%$ i1 = (d1==dates('1949Y'));
%$ i2 = (d2==dates('1950Q1'));
%$ i3 = (d3==dates('1949Q4'));
%$ i4 = (d4==dates('1950M1'));
%$ i5 = (d5==dates('1949M12'));
%$
%$ % Check the results.
%$ t(1) = dassert(i1,true);
%$ t(2) = dassert(i2,true);
%$ t(3) = dassert(i3,true);
%$ t(4) = dassert(i4,true);
%$ t(5) = dassert(i5,true);
%$ T = all(t);
%@eof:1

%@test:2
%$ d1 = dates('1950Q1','1950Q2','1950Q3','1950Q4','1951Q1');
%$ d2 = dates('1949Q4','1950Q1','1950Q2','1950Q3','1950Q4');
%$ try
%$   d3 = -d1;
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(all(d2==d3),true);
%$ end
%$
%$ T = all(t);
%@eof:2