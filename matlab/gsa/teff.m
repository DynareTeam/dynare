function [yt, j0, ir, ic]=teff(T,Nsam,istable)
% [yt, j0, ir, ic]=teff(T,Nsam,istable)
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

ndim = (length(size(T)));
if ndim==3,
if nargin==1,
  Nsam=size(T,3);
  istable = [1:Nsam]';
end
tmax=max(T,[],3);
tmin=min(T,[],3);
[ir, ic]=(find( (tmax-tmin)>1.e-8));
j0 = length(ir);
yt=zeros(Nsam, j0);

for j=1:j0,
  y0=squeeze(T(ir(j),ic(j),:));
  %y1=ones(size(lpmat,1),1)*NaN;
  y1=ones(Nsam,1)*NaN;
  y1(istable,1)=y0;
  yt(:,j)=y1;
end

else
tmax=max(T,[],2);
tmin=min(T,[],2);
ir=(find( (tmax-tmin)>1.e-8));
j0 = length(ir);
yt=NaN(Nsam, j0);
yt(istable,:)=T(ir,:)';

  
end
%clear y0 y1;
