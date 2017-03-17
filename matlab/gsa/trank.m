function yr = trank(y);
% yr = trank(y);
% yr is the rank transformation of y
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu
%
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.

% Copyright (C) 2012 European Commission
% Copyright (C) 2012 Dynare Team
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

[nr, nc] = size(y);
for j=1:nc,
  [dum, is]=sort(y(:,j));
  yr(is,j)=[1:nr]'./nr;
end
