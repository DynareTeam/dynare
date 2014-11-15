function b = isquaterly(str)  % --*-- Unitary tests --*--

% Tests if the input can be interpreted as a quaterly date.
%
% INPUTS 
%  o str     string.
%
% OUTPUTS 
%  o b       integer scalar, equal to 1 if str can be interpreted as a quaterly date or 0 otherwise.

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

if ischar(str)
    if isempty(regexp(str,'^-?[0-9]*[Qq][1-4]$','once'))
        b = 0;
    else
        b = 1;
    end
else
    b = 0;
end

%@test:1
%$
%$ date_1 = '1950Q2';
%$ date_2 = '1950q2';
%$ date_3 = '-1950q2';
%$ date_4 = '1950q12';
%$ date_5 = '1950 azd ';
%$ date_6 = '1950Y';
%$ date_7 = '1950m24';
%$
%$ t(1) = dassert(isquaterly(date_1),1);
%$ t(2) = dassert(isquaterly(date_2),1);
%$ t(3) = dassert(isquaterly(date_3),1);
%$ t(4) = dassert(isquaterly(date_4),0);
%$ t(5) = dassert(isquaterly(date_5),0);
%$ t(6) = dassert(isquaterly(date_6),0);
%$ t(7) = dassert(isquaterly(date_7),0);
%$ T = all(t);
%@eof:1