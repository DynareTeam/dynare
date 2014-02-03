function [ldens,Dldens,D2ldens] = lpdfgbeta(x,a,b,aa,bb);
% Evaluates the logged BETA PDF at x. 
%
% INPUTS 
%    x     [double]  m*n matrix of loactions,
%    a     [double]  m*n matrix of First BETA distribution parameters, 
%    b     [double]  m*n matrix of Second BETA distribution parameters, 
%    aa    [double]  m*n matrix of lower bounds for (generalized) distribution, 
%    bb    [double]  m*n matrix of upper bounds for (generalized) distribution
%
% OUTPUTS 
%    ldens [double]  m*n matrix of logged (generalized) BETA densities.
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2011 Dynare Team
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

ldens = -Inf( size(x) ) ;
idx = find( (x-aa)>0 & (x-bb)<0 ) ;

if length(a)==1
    ldens(idx) = -betaln(a,b) + (a-1)*log(x(idx)-aa) + (b-1)*log(bb-x(idx)) - (a+b-1)*log(bb-aa) ;
else
    ldens(idx) = -betaln(a(idx),b(idx)) + (a(idx)-1).*log(x(idx)-aa(idx)) + (b(idx)-1).*log(bb(idx)-x(idx)) - (a(idx)+b(idx)-1).*log(bb(idx)-aa(idx));
end


if nargout >1 
    if length(a)==1
        Dldens(idx) = (a-1)./(x(idx)-aa) - (b-1)./(bb-x(idx)) ;
    else
        Dldens(idx) = (a(idx)-1)./(x(idx)-aa(idx)) - (b(idx)-1)./(bb(idx)-x(idx));
    end
end


if nargout == 3 
    if length(a)==1
        D2ldens(idx) = -(a-1)./(x(idx)-aa).^2 - (b-1)./(bb-x(idx)).^2 ;
    else
        D2ldens(idx) = -(a(idx)-1)./(x(idx)-aa(idx)).^2 - (b(idx)-1)./(bb(idx)-x(idx)).^2;
    end
end