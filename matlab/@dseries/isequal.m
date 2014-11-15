function C = isequal(A, B, tol)

% Overloads the isequal Matlab/Octave's function.
%
% INPUTS 
%  o A      dseries object (T periods, N variables).
%  o B      dseries object (T periods, N variables).
%  o tol    tolerance parameter.
%
% OUTPUTS 
%  o C      Integer scalar equal to zero or one.

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

if nargin~=2
    error('dseries::isequal: I need exactly two input arguments!')
end

if ~isdseries(B)
    error('dseries::isequal: Both input arguments must be dseries objects!')
end

if ~isequal(nobs(A), nobs(B))
    C = 0;
    return
end

if ~isequal(vobs(A), vobs(B))
    C = 0;
    return
end

if ~isequal(frequency(A),frequency(B))
    C = 0;
    return
end

if ~isequal(A.dates,B.dates)
    C = 0;
    return
end

if ~isequal(A.name,B.name)
    warning('dseries::isequal: Both input arguments do not have the same variables!')
end

if ~isequal(A.tex,B.tex)
    warning('dseries::isequal: Both input arguments do not have the same tex names!')
end

if nargin<3
    C = isequal(A.data, B.data);
else
    C = ~(max(abs(A.data(:)-B.data(:)))>tol);
end