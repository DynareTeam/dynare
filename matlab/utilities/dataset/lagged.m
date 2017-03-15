function xlag = lagged(x, n);
% xlag = lagged(x, n);
% applies n-lags backward shift operator to x
% 
% INPUTS
% x    = time series 
% n    = number of backward shifts [DEFAULT=1]
% 
% OUTPUT
% xlag = backward shifted series

% Copyright (C) 2017 Dynare Team
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

if nargin==1, n=1; end

x=x(:);
xlag=[NaN(n,1); x(1:end-n)];
