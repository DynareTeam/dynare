function [Z,ST,R1,QT,Pstar,Pinf] = schur_statespace_transformation(mf,T,R,Q,qz_criterium, restrict_columns)
% function [Z,ST,QT,R1,Pstar,Pinf] = schur_statespace_transformation(mf,T,R,Q,qz_criterium, restrict_columns)
% Kitagawa transformation of state space system with a quasi-triangular
% transition matrix with unit roots at the top, but excluding zero columns of the transition matrix. 
% Computation of Pstar and Pinf for Durbin and Koopman Diffuse filter
%
% The transformed state space is 
%     y = [ss; z; x];
%     s = static variables (zero columns of T)
%     z = unit roots
%     x = stable roots
%     ss = s - z = stationarized static variables
% 
% INPUTS 
%   mf           [integer]    vector of indices of observed variables in
%                             state vector
%   T            [double]     matrix of transition
%   R            [double]     matrix of structural shock effects
%   Q            [double]     matrix of covariance of structural shocks
%   qz_criterium [double]     numerical criterium for unit roots   
%  
% OUTPUTS 
%   Z            [double]     transformed matrix of measurement equation
%   ST           [double]     tranformed matrix of transition
%   R1           [double]     tranformed matrix of structural shock effects
%   QT           [double]     matrix of Schur vectors
%   Pstar        [double]     matrix of covariance of stationary part
%   Pinf         [double]     matrix of covariance initialization for
%                             nonstationary part    
%    
% ALGORITHM 
%   Real Schur transformation of transition equation
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2006-2011 Dynare Team
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

np = size(T,1);
if nargin == 6,
    indx = restrict_columns;
    indx0=find(~ismember([1:np],indx));
else
    indx=(find(max(abs(T))>0));
    indx0=(find(max(abs(T))==0));
end
T0=T(indx0,indx); %static variables vs. dynamic ones
R0=R(indx0,:); % matrix of shocks for static variables

% perform Kitagawa transformation only for non-zero columns of T
T=T(indx,indx);
R=R(indx,:);
np = size(T,1);
[QT,ST] = schur(T);
e1 = abs(ordeig(ST)) > 2-qz_criterium;
[QT,ST] = ordschur(QT,ST,e1);
k = find(abs(ordeig(ST)) > 2-qz_criterium);
nk = length(k);
nk1 = nk+1;
Pstar = zeros(np,np);
B = QT'*R*Q*R'*QT;
i = np;
while i >= nk+2
    if ST(i,i-1) == 0
        if i == np
            c = zeros(np-nk,1);
        else
            c = ST(nk1:i,:)*(Pstar(:,i+1:end)*ST(i,i+1:end)')+...
                ST(i,i)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i);
        end
        q = eye(i-nk)-ST(nk1:i,nk1:i)*ST(i,i);
        Pstar(nk1:i,i) = q\(B(nk1:i,i)+c);
        Pstar(i,nk1:i-1) = Pstar(nk1:i-1,i)';
        i = i - 1;
    else
        if i == np
            c = zeros(np-nk,1);
            c1 = zeros(np-nk,1);
        else
            c = ST(nk1:i,:)*(Pstar(:,i+1:end)*ST(i,i+1:end)')+...
                ST(i,i)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i)+...
                ST(i,i-1)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i-1);
            c1 = ST(nk1:i,:)*(Pstar(:,i+1:end)*ST(i-1,i+1:end)')+...
                 ST(i-1,i-1)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i-1)+...
                 ST(i-1,i)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i);
        end
        q = [eye(i-nk)-ST(nk1:i,nk1:i)*ST(i,i) -ST(nk1:i,nk1:i)*ST(i,i-1);...
             -ST(nk1:i,nk1:i)*ST(i-1,i) eye(i-nk)-ST(nk1:i,nk1:i)*ST(i-1,i-1)];
        z =  q\[B(nk1:i,i)+c;B(nk1:i,i-1)+c1];
        Pstar(nk1:i,i) = z(1:(i-nk));
        Pstar(nk1:i,i-1) = z(i-nk+1:end);
        Pstar(i,nk1:i-1) = Pstar(nk1:i-1,i)';
        Pstar(i-1,nk1:i-2) = Pstar(nk1:i-2,i-1)';
        i = i - 2;
    end
end
if i == nk+1
    c = ST(nk+1,:)*(Pstar(:,nk+2:end)*ST(nk1,nk+2:end)')+ST(nk1,nk1)*ST(nk1,nk+2:end)*Pstar(nk+2:end,nk1);
    Pstar(nk1,nk1)=(B(nk1,nk1)+c)/(1-ST(nk1,nk1)*ST(nk1,nk1));
end

R1 = QT'*R;
ST1=ST;

% now I recover stationarized static variables
% using 
% ss = s-z and
% z-z(-1) (growth rates of unit roots) only depends on stationary variables
np0=length(indx0);
Pstar = blkdiag(zeros(np0),Pstar);
ST = [zeros(length(Pstar),length(indx0)) [T0*QT ;ST]];
R1 = [R0; R1];
ST0=ST;
ST0(:,1:np0+nk)=0;
ST0(np0+1:np0+nk,:)=0;
ST0(1:np0,np0+nk+1:end) = ST(1:np0,np0+nk+1:end)-ST(1:np0,np0+1:np0+nk)*ST(np0+1:np0+nk,np0+nk+1:end);
R10 = R1;
R10(np0+1:np0+nk,:)=0;
R10(1:np0,:) = R1(1:np0,:)-ST(1:np0,np0+1:np0+nk)*R1(np0+1:np0+nk,:);
Pstar = ...
     ST0*Pstar*ST0'+ ...
         R10*Q*R10';
QT = blkdiag(eye(np0),QT);
QT(1:np0,np0+1:np0+nk) = QT(1:np0,np0+1:np0+nk)+ST(1:np0,np0+1:np0+nk);
% reorder QT entries 
QT([indx0(:); indx(:)],:) = QT; 
Z = QT(mf,:);
ST(1:np0,:) = ST0(1:np0,:);
R1(1:np0,:) = R10(1:np0,:);

% stochastic trends with no influence on observed variables are
% arbitrarily initialized to zero
Pinf = zeros(np,np);
Pinf(1:nk,1:nk) = eye(nk);
[QQ,RR,EE] = qr(Z*ST(:,1+np0:nk+np0),0);
k = find(abs(diag([RR; zeros(nk-size(Z,1),size(RR,2))])) < 1e-8);
if length(k) > 0
    k1 = EE(:,k);
    dd =ones(nk,1);
    dd(k1) = zeros(length(k1),1);
    Pinf(1:nk,1:nk) = diag(dd);
end
Pinf = blkdiag(zeros(np0),Pinf);

