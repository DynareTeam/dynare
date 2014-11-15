function C = compare_vectors(f, A, B)  % --*-- Unitary tests --*--
 
% Performs lexicographical comparison of vectors.
%
% INPUTS
%  o f    function handle: @lt, @gt, @le or @ge
%  o A    vector of real numbers.
%  o B    vector of real numbers.
%
% OUTPUTS 
%  o C    integer scalar, 1 or 0.
%
% REMARKS 
%  o It is assumed that vectors A and B have the same number of elements.
    
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

if ~isvector(A) || ~isvector(B)
    error('compare_vectors:: Input arguments a and b must be vectors!')
end

if ~isequal(length(A),length(B))
    error('compare_vectors_lt:: Input arguments a and b must be of same length!')
end

fstr = func2str(f);

if ~ismember(fstr, {'lt', 'gt', 'le', 'ge'})
    error('compare_vectors:: First input argument must be one of the following function handles: @lt, @gt, @le or @ge!')
end

strict_inequality = ismember(fstr, {'gt','lt'});

if isequal(length(A),1)
    if feval(f, A, B)
        C = 1;
    else
        C = 0;
    end
else
    if feval(f, A(1), B(1))
        if strict_inequality
            C = 1;
        else
            if isequal(A(1),B(1))
                C = compare_vectors(f, A(2:end), B(2:end));
            else
                C = 1;
            end
        end
    else
        if strict_inequality
            if isequal(A(1),B(1))
                C = compare_vectors(f, A(2:end), B(2:end));
            else
                C = 0;
            end
        else
            C = 0;
        end
    end
end

%@test:1
%$ t(1) = dassert(compare_vectors(@lt, [1990 3], [1991 3]), 1);
%$ t(2) = dassert(compare_vectors(@lt, [1990 3], [1990 3]), 0);
%$ t(3) = dassert(compare_vectors(@le, [1990 3], [1990 3]), 1);
%$ t(4) = dassert(compare_vectors(@lt, [1990 3], [1990 4]), 1);
%$ t(5) = dassert(compare_vectors(@le, [1990 3], [1990 4]), 1);
%$ t(6) = dassert(compare_vectors(@gt, [1990 3], [1991 3]), 0);
%$ t(7) = dassert(compare_vectors(@gt, [1990 3], [1990 3]), 0);
%$ t(8) = dassert(compare_vectors(@ge, [1990 3], [1990 3]), 1);
%$ t(9) = dassert(compare_vectors(@gt, [1990 3], [1990 4]), 0);
%$ t(10) = dassert(compare_vectors(@ge, [1990 3], [1990 4]), 0);
%$ t(11) = dassert(compare_vectors(@le, [1991 3], [1990 4]), 0);
%$ t(12) = dassert(compare_vectors(@le, [1991 3], [1990 2]), 0);
%$ t(13) = dassert(compare_vectors(@le, [1945 2], [1950, 1]),1);
%$ T = all(t);
%@eof:1