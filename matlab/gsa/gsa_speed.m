function [tadj, iff] = gsa_speed(A,B,mf,p),
% [tadj, iff] = gsa_speed(A,B,mf,p),
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

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

nvar=length(mf);
nstate= size(A,1);
nshock = size(B,2);
nrun=size(B,3);

iff=zeros(nvar,nshock,nrun);
tadj=iff;
disp('Computing speed of adjustement ...')
h = dyn_waitbar(0,'Speed of adjustement...');

for i=1:nrun,
  irf=zeros(nvar,nshock);
  a=squeeze(A(:,:,i));
  b=squeeze(B(:,:,i));
  IFF=inv(eye(nstate)-a)*b;
  iff(:,:,i)=IFF(mf,:);
  IF=IFF-b;
  
  t=0;
  while any(any(irf<0.5))
    t=t+1;
    IFT=((eye(nstate)-a^(t+1))*inv(eye(nstate)-a))*b-b;
    irf=IFT(mf,:)./(IF(mf,:)+eps);
    irf = irf.*(abs(IF(mf,:))>1.e-7)+(abs(IF(mf,:))<=1.e-7);
    %irf=ft(mf,:);
    tt=(irf>0.5).*t;
    tadj(:,:,i)=((tt-tadj(:,:,i))==tt).*tt+tadj(:,:,i);
  end
  dyn_waitbar(i/nrun,h)
end
skipline()
disp('.. done !')
dyn_waitbar_close(h)
