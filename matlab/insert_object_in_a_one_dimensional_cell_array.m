function one_dimensional_cell_array = insert_object_in_a_one_dimensional_cell_array(one_dimensional_cell_array, object, i) % --*-- Unitary tests --*--

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

if nargin<2
    error('insert_object_in_a_one_dimensional_cell_array: I need at least two input arguments!')
end

if ~iscell(one_dimensional_cell_array)
    error(['insert_column_vector_in_a_matrix: First input ''' inputname(1) ''' must be a cell array!'])
end

[nr, nc] = size(one_dimensional_cell_array);

if ~isequal(nr*nc,numel(one_dimensional_cell_array))
    error(['insert_column_vector_in_a_matrix: First input ''' inputname(1) ''' must be a one dimensional cell array!'])
end

n = numel(one_dimensional_cell_array);

if nargin<3
    i = n+1;
end

one_dimensional_cell_array = one_dimensional_cell_array(:);

switch i
  case n+1
    one_dimensional_cell_array = [one_dimensional_cell_array; object];
  case 1
    one_dimensional_cell_array = [object; one_dimensional_cell_array];
  otherwise
    one_dimensional_cell_array = [one_dimensional_cell_array(1:i-1); object; one_dimensional_cell_array(i:n)];
end

if nc>nr
    one_dimensional_cell_array = transpose(one_dimensional_cell_array);
end

%@test:1
%$ A = {'A1'; 'A2'; 'A3'}; b = [1, pi];
%$
%$ try
%$   C4 = insert_object_in_a_one_dimensional_cell_array(A, b);
%$   C1 = insert_object_in_a_one_dimensional_cell_array(A, b, 1);
%$   C2 = insert_object_in_a_one_dimensional_cell_array(A, b, 2);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(C4,[A;b]);
%$   t(3) = dassert(C1,[b;A],1e-15);
%$   t(4) = dassert(C2,[A(1); b; A(2:3)]);
%$ end
%$ T = all(t);
%@eof:1