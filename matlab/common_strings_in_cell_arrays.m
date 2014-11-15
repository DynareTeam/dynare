function [n,message] = common_strings_in_cell_arrays(a,b)

%@info:
%! @deftypefn {Function File} {[@var{n}, @var{message} ] =} common_strings_in_cell_arrays (@var{a}, @var{b})
%! @anchor{common_strings_in_cell_arrays}
%! @sp 1
%! Counts the number of common strings in two cell arrays of strings @var{a} and @var{b}. 
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! One dimensional cell array of strings.
%! @item b
%! One dimensional cell array of strings.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item n
%! Scalar integer, number of common strings.
%! @item message
%! String, formatted list of common strings. 
%! @end table
%! @end deftypefn
%@eod:

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

common_strings = intersect(a,b);
n = numel(common_strings);

if nargout>1
    switch n
      case 0
        message = [];
      case 1
        message = common_strings{1};
      case 2
        message = [common_strings{1} ' and ' common_strings{2}];
      otherwise
        message = common_strings{1};
        for i = 2:n-1
            message = [message ', ' common_strings{i}];
        end
        message = [message ' and ' common_strings{n}];
    end
end

%@test:1
%$ % Define cell arrays of strings.
%$ A = {'A1'; 'A2'; 'A3'; 'A4'};
%$ B = {'B1'; 'B2'; 'B3'; 'B4'};
%$
%$ try
%$    n = common_strings_in_cell_arrays(A,B);
%$    [m,message] = common_strings_in_cell_arrays(A,B);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(n,m);
%$    t(3) = dassert(n,0);
%$    t(4) = dassert(isempty(message),1);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define cell arrays of strings.
%$ A = {'A1'; 'A2'; 'A3'; 'A4'};
%$ B = {'B1'; 'A2'; 'B3'; 'A4'};
%$
%$ try
%$    n = common_strings_in_cell_arrays(A,B);
%$    [m,message] = common_strings_in_cell_arrays(A,B);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(n,m);
%$    t(3) = dassert(n,2);
%$    t(4) = dassert(message,'A2 and A4');
%$ end
%$ T = all(t);
%@eof:2

%@test:2
%$ % Define cell arrays of strings.
%$ A = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'A6'};
%$ B = {'B1'; 'A2'; 'B3'; 'A4'; 'A1'};
%$
%$ try
%$    n = common_strings_in_cell_arrays(A,B);
%$    [m,message] = common_strings_in_cell_arrays(A,B);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dassert(n,m);
%$    t(3) = dassert(n,3);
%$    t(4) = dassert(message,'A2, A4 and A1');
%$ end
%$ T = all(t);
%@eof:3

