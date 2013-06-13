function c = plus(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} plus (@var{a},@var{b})
%! @anchor{@dynDate/plus}
%! @sp 1
%! Overloads the plus (addition) operator for the Dynare dates class (@ref{dynDate}). Given an initial date @var{a},
%! computes a new date @var{c} by adding the number of periods @var{b}.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Dynare date object instantiated by @ref{dynDate}.
%! @item b
%! Positive scalar integer, the number of periods
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item c
%! Dynare date object instantiated by @ref{dynDate}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{@@dynDate/eq}
%!
%! @end deftypefn
%@eod:

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

if ~isa(a,'dynDate')
    error(['dynDate::plus: Input argument ' inputname(1) ' must be a dynDate object!'])
end

if b<0 || ~isint(b)
    error(['dynDate::plus: Input argument ' inputname(2) ' must be a positive integer'])
end

if b==0
    c = a;
    return
end

switch a.freq
  case 1
    c = a;
    c.time(1) = a.time(1) + b - 1;
  case {4,12,52}
    c = a;
    n1 = b + a.time(2);
    n2 = floor((n1 - 1)/a.freq);
    n3 = mod(n1 - 1,a.freq) + 1;
    c.time(2) = n3;
    c.time(1) = c.time(1)+n2;
  otherwise
    error('dynDate::plus: Unknown frequency!')
end

%@test:1
%$ % Define some dates
%$ date_1 = 1950;
%$ date_2 = '1950Q4';
%$ date_3 = '2000M3';
%$
%$ % Call the tested routine.
%$ d_1 = dynDate(date_1);
%$ d_2 = dynDate(date_2);
%$ d_3 = dynDate(date_3);
%$
%$ d1 = d_1+3;
%$ d2 = d_2+5;
%$ d3 = d_3+15;
%$ d4 = d_3+10;
%$
%$ % Expected results.
%$ e1 = dynDate(1952);
%$ e2 = dynDate('1951Q4');
%$ e3 = dynDate('2001M5');
%$ e4 = dynDate('2000M12');
%$
%$ % Check the results.
%$ t(1) = dyn_assert(e1.time,d1.time);
%$ t(2) = dyn_assert(e2.time,d2.time);
%$ t(3) = dyn_assert(e3.time,d3.time);
%$ t(4) = dyn_assert(e4.time,d4.time);
%$ T = all(t);
%@eof:1
