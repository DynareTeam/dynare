function [y, meany, stdy] = stand_(x)
% STAND_  Standardise a matrix by columns
%
% [x,my,sy]=stand(y)
%
% y: Time series (column matrix)
%
% x: standardised equivalent of y
% my: Vector of mean values for each column of y
% sy: Vector of standard deviations for each column of y
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2012 European Commission
% Copyright (C) 2012 Dynare Team%
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

if nargin==0,
    return;
end

for j=1:size(x,2);
meany(j)=mean(x(find(~isnan(x(:,j))),j));
stdy(j)=std(x(find(~isnan(x(:,j))),j));
    y(:,j)=(x(:,j)-meany(j))./stdy(j);
end
% end of m-file