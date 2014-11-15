function b = isdate(str)  % --*-- Unitary tests --*--

% Tests if the input string can be interpreted as a date.
%
% INPUTS 
%  o str     string.
%
% OUTPUTS 
%  o b       integer scalar, equal to 1 if str can be interpreted as a date or 0 otherwise.

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

if isnumeric(str) && isscalar(str)
    b = true;
    return
end

b = isstringdate(str);

%@test:1
%$
%$ date_1 = 1950;
%$ date_2 = '1950m2';
%$ date_3 = '-1950m2';
%$ date_4 = '1950m52';
%$ date_5 = ' 1950';
%$ date_6 = '1950Y';
%$ date_7 = '-1950a';
%$ date_8 = '1950m ';
%$
%$ t(1) = dassert(isdate(date_1),true);
%$ t(2) = dassert(isdate(date_2),true);
%$ t(3) = dassert(isdate(date_3),true);
%$ t(4) = dassert(isdate(date_4),false);
%$ t(5) = dassert(isdate(date_5),false);
%$ t(6) = dassert(isdate(date_6),true);
%$ t(7) = dassert(isdate(date_7),true);
%$ t(8) = dassert(isdate(date_8),false);
%$ T = all(t);
%@eof:1