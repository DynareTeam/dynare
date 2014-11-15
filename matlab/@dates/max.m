function C = max(varargin)
    
% Overloads the max function for dates objects.
    
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

switch nargin
  case 1
    switch length(varargin{1})
      case 0
        C= dates();
      case 1
        C = varargin{1};
      otherwise
        tmp = sortrows(varargin{1}.time);
        C = dates();
        C.freq = varargin{1}.freq;
        C.ndat = 1;
        C.time = tmp(varargin{1}.ndat,:);
    end
  otherwise
    C = max(horzcat(varargin{:}));
end

%@test:1
%$ % Define some dates
%$ d3 = dates('1950q2');
%$ d4 = dates('1950Q3');
%$ d5 = dates('1950m1');
%$ d6 = dates('1948M6');
%$ m2 = max(d3,d4);
%$ i2 = (m2==d4);
%$ m3 = max(d5,d6);
%$ i3 = (m3==d5);
%$
%$ % Check the results.
%$ t(1) = dassert(i2,1);
%$ t(2) = dassert(i3,1);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ d = dates('1950Q2','1951Q3','1949Q1','1950Q4');
%$ m = max(d);
%$ i = (m==dates('1951Q3'));
%$
%$ % Check the results.
%$ t(1) = dassert(i,1);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dates
%$ m = max(dates('1950Q2','1951Q3'),dates('1949Q1'),dates('1950Q4'));
%$ i = (m==dates('1951Q3'));
%$
%$ % Check the results.
%$ t(1) = dassert(i,1);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define some dates
%$ m = max(dates('1950Q2'),dates('1951Q3'),dates('1949Q1'),dates('1950Q4'));
%$ i = (m==dates('1951Q3'));
%$
%$ % Check the results.
%$ t(1) = dassert(i,1);
%$ T = all(t);
%@eof:4