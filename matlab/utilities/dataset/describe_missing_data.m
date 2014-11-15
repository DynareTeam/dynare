function [i,n,s,j] = describe_missing_data(data)
% This function reads the dataset and determines the location of the missing observations (defined by NaNs)

%@info:
%! @deftypefn {Function File} {[@var{i}, @var{n}, @var{s}, @var{j} ] =} describe_missing_data (@var{data}, @var{gend}, @var{nvarobs})
%! This function reads the dataset and determines where are the missing observations.
%!
%! @strong{Inputs}
%! @table @var
%! @item data
%! Real matrix (T-by-N) for the dataset.
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item i
%! cell array (1-by-T). Each element is a @math{p_t\times 1} column vector of indices targeting the non-NaN variables at time t.
%! @item n
%! Integer scalar. The effective number of observations:
%!    @math(n=\sum_{t=1}^T p_t)
%! @item s
%! Integer scalar. The value of the time index such that @math(p_t=p_s) for all @math(t\geq s).
%! @item j
%! cell array (1-by-N). Each element is a column vector targeting to the non-NaN observations of a variable.
%! @end table
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2008-2014 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

[observation_index,variable_index] = find(~isnan(data));
[T,N] = size(data);

i = cell(1,T);
j = cell(1,N);
missing_observations_counter = NaN(T,1);

for obs=1:T
    idx = find(observation_index==obs);
    tmp = variable_index(idx);
    missing_observations_counter(obs,1) = N-length(tmp);
    if rows(tmp(:))
        i(obs) = { tmp(:) };
    else
        i(obs) = { [] };
    end
end

missing_observations_counter = cumsum(missing_observations_counter);

n = length(variable_index);

if ~missing_observations_counter
    s = 1;
else
    tmp = find(missing_observations_counter>=(T*N-n));
    s = tmp(1)+1;
end

if nargout>3
    for var=1:N
        idx = find(variable_index==var);
        tmp = observation_index(idx);
        j(var) = { tmp(:) };
    end
end


%@test:1
%$ % Define a data set.
%$ A = [ 1    1   ;   ...
%$       1    NaN ;   ...
%$       NaN  1   ;   ...
%$       1    1   ;   ...
%$       NaN  NaN ;   ...
%$       1    NaN ;   ...
%$       1    NaN ;   ...
%$       1    1   ;   ...
%$       1    1   ;   ...
%$       1    1   ;   ...
%$       1    1  ];
%$
%$ % Define expected results.
%$ eB = cell(1,11);
%$ eB(1)  = { transpose(1:2) };
%$ eB(2)  = { 1 };
%$ eB(3)  = { 2 };
%$ eB(4)  = { transpose(1:2)};
%$ eB(5)  = { [] };
%$ eB(6)  = { 1 };
%$ eB(7)  = { 1 };
%$ eB(8)  = { transpose(1:2) };
%$ eB(9)  = { transpose(1:2) };
%$ eB(10) = { transpose(1:2) };
%$ eB(11) = { transpose(1:2) };
%$ eC = 16;
%$ eD = 8;
%$ eE = cell(1,2);
%$ eE(1) = { [1; 2; 4; transpose(6:11)] };
%$ eE(2) = { [1; 3; 4; transpose(8:11)] };
%$
%$ % Call the tested routine.
%$ [B,C,D,E] = describe_missing_data(transpose(A));
%$
%$ % Check the results.
%$ t(1) = dassert(B,eB);
%$ t(2) = dassert(C,eC);
%$ t(3) = dassert(D,eD);
%$ t(4) = dassert(E,eE);
%$ T = all(t);
%@eof:1