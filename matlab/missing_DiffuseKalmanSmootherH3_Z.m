function [alphahat,epsilonhat,etahat,a,P,aK,PK,decomp] = missing_DiffuseKalmanSmootherH3_Z(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,diffuse_kalman_tol,decomp_flag)
% function [alphahat,epsilonhat,etahat,a1,P,aK,PK,d,decomp] = missing_DiffuseKalmanSmootherH3_Z(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,decomp_flag)
% Computes the diffuse kalman smoother in the case of a singular var-cov matrix.
% Univariate treatment of multivariate time series.
%
% INPUTS
%    T:        mm*mm matrix     state transition matrix
%    Z:        pp*mm matrix     selector matrix for observables in augmented state vector
%    R:        mm*rr matrix     second matrix of the state equation relating the structural innovations to the state variables
%    Q:        rr*rr matrix     covariance matrix of structural errors
%    H:        pp*1             vector of variance of measurement errors
%    Pinf1:    mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:   mm*mm variance-covariance matrix with stationary variables
%    Y:        pp*1 vector
%    pp:       number of observed variables
%    mm:       number of state variables
%    smpl:     sample size
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    nk        number of forecasting periods
%    kalman_tol   tolerance for zero divider
%    diffuse_kalman_tol   tolerance for zero divider
%    decomp_flag  if true, compute filter decomposition
%
% OUTPUTS
%    alphahat: smoothed state variables (a_{t|T})
%    epsilonhat: measurement errors
%    etahat:   smoothed shocks
%    a:        matrix of updated variables (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t})
%              (meaningless for periods 1:d)
%    P:        3D array of one-step ahead forecast error variance
%              matrices
%    PK:       4D array of k-step ahead forecast error variance
%              matrices (meaningless for periods 1:d)
%    decomp:   decomposition of the effect of shocks on filtered values
%
% Notes:
%   Outputs are stored in decision-rule order, i.e. to get variables in order of declaration
%   as in M_.endo_names, ones needs code along the lines of:
%   variables_declaration_order(dr.order_var,:) = alphahat
% 
% SPECIAL REQUIREMENTS
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003), in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98. 

% Copyright (C) 2004-2016 Dynare Team
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

% Modified by M. Ratto
% New output argument aK: 1-step to nk-stpe ahed predictions)
% New input argument nk: max order of predictions in aK


if size(H,2)>1 
    error('missing_DiffuseKalmanSmootherH3_Z:: H is not a vector. This must not happens')
end
    
d = 0;
decomp = [];
spinf           = size(Pinf1);
spstar          = size(Pstar1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl);
a1              = zeros(mm,smpl+1);
aK              = zeros(nk,mm,smpl+nk);

Fstar           = zeros(pp,smpl);
Finf            = zeros(pp,smpl);
Fi              = zeros(pp,smpl);
Ki              = zeros(mm,pp,smpl);
Kstar           = zeros(mm,pp,smpl);
P               = zeros(mm,mm,smpl+1);
P1              = P;
PK              = zeros(nk,mm,mm,smpl+nk);
Pstar           = zeros(spstar(1),spstar(2),smpl); Pstar(:,:,1) = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl); Pinf(:,:,1) = Pinf1;
Pstar1          = Pstar;
Pinf1           = Pinf;
steady          = smpl;
rr              = size(Q,1); % number of structural shocks
QQ              = R*Q*transpose(R);
QRt                     = Q*transpose(R);
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
epsilonhat      = zeros(rr,smpl);
r               = zeros(mm,smpl);

t = 0;
icc=0;
newRank = rank(Pinf(:,:,1),diffuse_kalman_tol);
while newRank && t < smpl
    t = t+1;
    a(:,t) = a1(:,t);
    Pstar1(:,:,t) = Pstar(:,:,t);
    Pinf1(:,:,t) = Pinf(:,:,t);
    di = data_index{t}';
    for i=di
        Zi = Z(i,:);
        v(i,t)      = Y(i,t)-Zi*a(:,t);
        Fstar(i,t)  = Zi*Pstar(:,:,t)*Zi' +H(i);
        Finf(i,t)   = Zi*Pinf(:,:,t)*Zi';
        Kstar(:,i,t) = Pstar(:,:,t)*Zi';
        if Finf(i,t) > diffuse_kalman_tol && newRank
            icc=icc+1;
            Kinf(:,i,t)       = Pinf(:,:,t)*Zi';
            Kinf_Finf         = Kinf(:,i,t)/Finf(i,t);
            a(:,t)            = a(:,t) + Kinf_Finf*v(i,t);
            Pstar(:,:,t)      = Pstar(:,:,t) + ...
                Kinf(:,i,t)*Kinf_Finf'*(Fstar(i,t)/Finf(i,t)) - ...
                Kstar(:,i,t)*Kinf_Finf' - ...
                Kinf_Finf*Kstar(:,i,t)';
            Pinf(:,:,t)       = Pinf(:,:,t) - Kinf(:,i,t)*Kinf(:,i,t)'/Finf(i,t);
        elseif Fstar(i,t) > kalman_tol 
            a(:,t)            = a(:,t) + Kstar(:,i,t)*v(i,t)/Fstar(i,t);
            Pstar(:,:,t)      = Pstar(:,:,t) - Kstar(:,i,t)*Kstar(:,i,t)'/Fstar(i,t);
        end
    end
    if newRank
        oldRank = rank(Pinf(:,:,t),diffuse_kalman_tol);
    else
        oldRank = 0;
    end
    a1(:,t+1) = T*a(:,t);
    aK(1,:,t+1) = a1(:,t+1); 
    for jnk=2:nk
        aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
    end
    Pstar(:,:,t+1) = T*Pstar(:,:,t)*T'+ QQ;
    Pinf(:,:,t+1) = T*Pinf(:,:,t)*T';
    P0=Pinf(:,:,t+1);
    if newRank,
        newRank       = rank(Pinf(:,:,t+1),diffuse_kalman_tol);
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')   
    end
end


d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
Fstar = Fstar(:,1:d);
Finf = Finf(:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
Pstar1 = Pstar1(:,:,1:d);
Pinf1  = Pinf1(:,:,1:d);
notsteady = 1;
while notsteady && t<smpl
    t = t+1;
    a(:,t) = a1(:,t);
    P1(:,:,t) = P(:,:,t);
    di = data_index{t}';
    for i=di
        Zi = Z(i,:);
        v(i,t)  = Y(i,t) - Zi*a(:,t);
        Fi(i,t) = Zi*P(:,:,t)*Zi' + H(i);
        Ki(:,i,t) = P(:,:,t)*Zi';
        if Fi(i,t) > kalman_tol
            a(:,t) = a(:,t) + Ki(:,i,t)*v(i,t)/Fi(i,t);
            P(:,:,t) = P(:,:,t) - Ki(:,i,t)*Ki(:,i,t)'/Fi(i,t);
        end
    end
    a1(:,t+1) = T*a(:,t);
    Pf          = P(:,:,t);
    aK(1,:,t+1) = a1(:,t+1); 
    for jnk=1:nk
        Pf = T*Pf*T' + QQ;
        PK(jnk,:,:,t+jnk) = Pf;
        if jnk>1
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));    
        end
    end
    P(:,:,t+1) = T*P(:,:,t)*T' + QQ;
    %  notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<kalman_tol);
end
% $$$ P_s=tril(P(:,:,t))+tril(P(:,:,t),-1)';
% $$$ P1_s=tril(P1(:,:,t))+tril(P1(:,:,t),-1)';
% $$$ Fi_s = Fi(:,t);
% $$$ Ki_s = Ki(:,:,t);
% $$$ L_s  =Li(:,:,:,t);
% $$$ if t<smpl
% $$$   P  = cat(3,P(:,:,1:t),repmat(P_s,[1 1 smpl-t]));
% $$$   P1  = cat(3,P1(:,:,1:t),repmat(P1_s,[1 1 smpl-t]));
% $$$   Fi = cat(2,Fi(:,1:t),repmat(Fi_s,[1 1 smpl-t]));
% $$$   Li  = cat(4,Li(:,:,:,1:t),repmat(L_s,[1 1 smpl-t]));
% $$$   Ki  = cat(3,Ki(:,:,1:t),repmat(Ki_s,[1 1 smpl-t]));
% $$$ end
% $$$ while t<smpl
% $$$   t=t+1;
% $$$   a(:,t) = a1(:,t);
% $$$   di = data_index{t}';
% $$$   for i=di
% $$$     Zi = Z(i,:);
% $$$     v(i,t)      = Y(i,t) - Zi*a(:,t);
% $$$     if Fi_s(i) > kalman_tol
% $$$       a(:,t) = a(:,t) + Ki_s(:,i)*v(i,t)/Fi_s(i);
% $$$     end
% $$$   end
% $$$   a1(:,t+1) = T*a(:,t);
% $$$   Pf          = P(:,:,t);
% $$$   for jnk=1:nk,
% $$$       Pf = T*Pf*T' + QQ;
% $$$       aK(jnk,:,t+jnk) = T^jnk*a(:,t);
% $$$       PK(jnk,:,:,t+jnk) = Pf;
% $$$   end
% $$$ end
ri=zeros(mm,1);
t = smpl+1;
while t > d+1
    t = t-1;
    di = flipud(data_index{t})';
    for i = di
        if Fi(i,t) > kalman_tol
            ri = Z(i,:)'/Fi(i,t)*v(i,t)+ri-Ki(:,i,t)'*ri/Fi(i,t)*Z(i,:)';
        end
    end
    r(:,t) = ri;
    alphahat(:,t) = a1(:,t) + P1(:,:,t)*r(:,t);
    etahat(:,t) = QRt*r(:,t);
    ri = T'*ri;
end
if d
    r0 = zeros(mm,d); 
    r0(:,d) = ri;
    r1 = zeros(mm,d);
    for t = d:-1:1
        di = flipud(data_index{t})';
        for i = di
            if Finf(i,t) > diffuse_kalman_tol 
                r1(:,t) = Z(i,:)'*v(i,t)/Finf(i,t) + ...
                          (Kinf(:,i,t)'*Fstar(i,t)/Finf(i,t)-Kstar(:,i,t)')*r0(:,t)/Finf(i,t)*Z(i,:)' + ...
                          r1(:,t)-Kinf(:,i,t)'*r1(:,t)/Finf(i,t)*Z(i,:)';
                r0(:,t) = r0(:,t)-Kinf(:,i,t)'*r0(:,t)/Finf(i,t)*Z(i,:)';
            elseif Fstar(i,t) > kalman_tol % step needed whe Finf == 0
                r0(:,t) = Z(i,:)'/Fstar(i,t)*v(i,t)+r0(:,t)-(Kstar(:,i,t)'*r0(:,t))/Fstar(i,t)*Z(i,:)';
            end
        end
        alphahat(:,t) = a1(:,t) + Pstar1(:,:,t)*r0(:,t) + Pinf1(:,:,t)*r1(:,t);
        r(:,t)        = r0(:,t);
        etahat(:,t)   = QRt*r(:,t);
        if t > 1
            r0(:,t-1) = T'*r0(:,t);
            r1(:,t-1) = T'*r1(:,t);
        end
    end
end

if decomp_flag
    decomp = zeros(nk,mm,rr,smpl+nk);
    ZRQinv = inv(Z*QQ*Z');
    for t = max(d,1):smpl
        ri_d = zeros(mm,1);
        di = flipud(data_index{t})';
        for i = di
            if Fi(i,t) > kalman_tol
                ri_d = Z(i,:)'/Fi(i,t)*v(i,t)+ri_d-Ki(:,i,t)'*ri_d/Fi(i,t)*Z(i,:)';
            end
        end
        
        % calculate eta_tm1t
        eta_tm1t = QRt*ri_d;
        % calculate decomposition
        Ttok = eye(mm,mm); 
        AAA = P1(:,:,t)*Z'*ZRQinv*Z*R;
        for h = 1:nk
            BBB = Ttok*AAA;
            for j=1:rr
                decomp(h,:,j,t+h) = eta_tm1t(j)*BBB(:,j);
            end
            Ttok = T*Ttok;
        end
    end
end

epsilonhat = Y - Z*alphahat;
