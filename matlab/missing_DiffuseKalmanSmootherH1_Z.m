function [alphahat,epsilonhat,etahat,atilde,P,aK,PK,decomp,V] = missing_DiffuseKalmanSmootherH1_Z(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,diffuse_kalman_tol,decomp_flag,state_uncertainty_flag)

% function [alphahat,epsilonhat,etahat,a,aK,PK,decomp] = DiffuseKalmanSmoother1(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,diffuse_kalman_tol,decomp_flag,state_uncertainty_flag)
% Computes the diffuse kalman smoother without measurement error, in the case of a non-singular var-cov matrix.
%
% INPUTS
%    T:        mm*mm matrix
%    Z:        pp*mm matrix
%    R:        mm*rr matrix
%    Q:        rr*rr matrix
%    H:        pp*pp matrix variance of measurement errors
%    Pinf1:    mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:   mm*mm variance-covariance matrix with stationary variables
%    Y:        pp*1 vector
%    pp:       number of observed variables
%    mm:       number of state variables
%    smpl:     sample size
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    nk        number of forecasting periods
%    kalman_tol   tolerance for reciprocal condition number
%    diffuse_kalman_tol   tolerance for reciprocal condition number (for Finf) and the rank of Pinf
%    decomp_flag  if true, compute filter decomposition
%    state_uncertainty_flag     if true, compute uncertainty about smoothed
%                               state estimate
%
% OUTPUTS
%    alphahat: smoothed variables (a_{t|T})
%    epsilonhat:smoothed measurement errors
%    etahat:   smoothed shocks
%    atilde:   matrix of updated variables (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t)
%              (meaningless for periods 1:d)
%    P:        3D array of one-step ahead forecast error variance
%              matrices
%    PK:       4D array of k-step ahead forecast error variance
%              matrices (meaningless for periods 1:d)
%    decomp:   decomposition of the effect of shocks on filtered values
%    V:        3D array of state uncertainty matrices
%
% Notes:
%   Outputs are stored in decision-rule order, i.e. to get variables in order of declaration
%   as in M_.endo_names, ones needs code along the lines of:
%   variables_declaration_order(dr.order_var,:) = alphahat
%
% SPECIAL REQUIREMENTS
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%   Durbin/Koopman (2012): "Time Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 5

% Copyright (C) 2004-2018 Dynare Team
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

% modified by M. Ratto:
% new output argument aK (1-step to k-step predictions)
% new options_.nk: the max step ahed prediction in aK (default is 4)
% new crit1 value for rank of Pinf
% it is assured that P is symmetric

d = 0;
decomp = [];
spinf           = size(Pinf1);
spstar          = size(Pstar1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl+1);
atilde          = zeros(mm,smpl);
aK              = zeros(nk,mm,smpl+nk);
PK              = zeros(nk,mm,mm,smpl+nk);
iF              = zeros(pp,pp,smpl);
Fstar           = zeros(pp,pp,smpl);
iFstar          = zeros(pp,pp,smpl);
iFinf           = zeros(pp,pp,smpl);
K               = zeros(mm,pp,smpl);
L               = zeros(mm,mm,smpl);
Linf            = zeros(mm,mm,smpl);
Lstar           = zeros(mm,mm,smpl);
Kstar           = zeros(mm,pp,smpl);
Kinf            = zeros(mm,pp,smpl);
P               = zeros(mm,mm,smpl+1);
Pstar           = zeros(spstar(1),spstar(2),smpl+1);
Pstar(:,:,1)    = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl+1);
Pinf(:,:,1)     = Pinf1;
rr              = size(Q,1);
QQ              = R*Q*transpose(R);
QRt             = Q*transpose(R);
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
epsilonhat      = zeros(rr,smpl);
r               = zeros(mm,smpl+1);
Finf_singular   = zeros(1,smpl);
if state_uncertainty_flag
    V               = zeros(mm,mm,smpl);
    N               = zeros(mm,mm,smpl+1);
else
    V=[];
end

t = 0;
while rank(Pinf(:,:,t+1),diffuse_kalman_tol) && t<smpl
    t = t+1;
    di = data_index{t};
    if isempty(di)
        %no observations, propagate forward without updating based on
        %observations
        atilde(:,t)     = a(:,t);
        a(:,t+1)        = T*atilde(:,t);
        Linf(:,:,t)     = T;
        Pstar(:,:,t+1)  = T*Pstar(:,:,t)*T' + QQ;
        Pinf(:,:,t+1)   = T*Pinf(:,:,t)*T';
    else
        ZZ = Z(di,:);                                                       %span selector matrix
        v(di,t)= Y(di,t) - ZZ*a(:,t);                                       %get prediction error v^(0) in (5.13) DK (2012)
        Finf = ZZ*Pinf(:,:,t)*ZZ';                                          % (5.7) in DK (2012)
        if rcond(Finf) < diffuse_kalman_tol                                 %F_{\infty,t} = 0
            if ~all(abs(Finf(:)) < diffuse_kalman_tol)                      %rank-deficient but not rank 0
                                                                            % The univariate diffuse kalman filter should be used.
                alphahat = Inf;
                return
            else                                                            %rank of F_{\infty,t} is 0
                Finf_singular(1,t) = 1;
                Fstar(di,di,t)  = ZZ*Pstar(:,:,t)*ZZ' + H(di,di);             % (5.7) in DK (2012)
                if rcond(Fstar(di,di,t)) < kalman_tol                         %F_{*} is singular
                    if ~all(all(abs(Fstar(di,di,t))<kalman_tol))
                        % The univariate diffuse kalman filter should be used.
                        alphahat = Inf;
                        return
                    else %rank 0
                        a(:,t+1) = T*a(:,t);
                        Pstar(:,:,t+1) = T*Pstar(:,:,t)*transpose(T)+QQ;
                        Pinf(:,:,t+1)  = T*Pinf(:,:,t)*transpose(T);
                    end
                else
                    iFstar(di,di,t) = inv(Fstar(di,di,t));
                    Kstar(:,di,t)   = Pstar(:,:,t)*ZZ'*iFstar(di,di,t);       %(5.15) of DK (2012) with Kstar=T^{-1}*K^(0)
                    Pinf(:,:,t+1)   = T*Pinf(:,:,t)*transpose(T);           % DK (2012), 5.16
                    Lstar(:,:,t)    = T - T*Kstar(:,di,t)*ZZ;               %L^(0) in DK (2012), eq. 5.12
                    Pstar(:,:,t+1)  = T*Pstar(:,:,t)*Lstar(:,:,t)'+QQ;      % (5.17) DK (2012)
                    a(:,t+1)        = T*(a(:,t)+Kstar(:,di,t)*v(di,t));       % (5.13) DK (2012)
                end
            end
        else
            %see notes in kalman_filter_d.m for details of computations
            iFinf(di,di,t)  = inv(Finf);
            Kinf(:,di,t)    = Pinf(:,:,t)*ZZ'*iFinf(di,di,t);               %define Kinf=T^{-1}*K_0 with M_{\infty}=Pinf*Z'
            atilde(:,t)     = a(:,t) + Kinf(:,di,t)*v(di,t);
            Linf(:,:,t)     = T - T*Kinf(:,di,t)*ZZ;                        %L^(0) in DK (2012), eq. 5.12
            Fstar(di,di,t)  = ZZ*Pstar(:,:,t)*ZZ' + H(di,di);               %(5.7) DK(2012)
            Kstar(:,di,t)   = (Pstar(:,:,t)*ZZ'-Kinf(:,di,t)*Fstar(di,di,t))*iFinf(di,di,t); %(5.12) DK(2012) with Kstar=T^{-1}*K^(1); note that there is a typo in DK (2003) with "+ Kinf" instead of "- Kinf", but it is correct in their appendix
            Pstar(:,:,t+1)  = T*Pstar(:,:,t)*Linf(:,:,t)'-T*Kinf(:,di,t)*Finf*Kstar(:,di,t)'*T' + QQ; %(5.14) DK(2012)
            Pinf(:,:,t+1)   = T*Pinf(:,:,t)*Linf(:,:,t)';                   %(5.14) DK(2012)
        end
        a(:,t+1)            = T*atilde(:,t);
        aK(1,:,t+1)         = a(:,t+1);
        % isn't a meaningless as long as we are in the diffuse part? MJ
        for jnk=2:nk
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        end
    end
end
d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
iFinf = iFinf(:,:,1:d);
iFstar= iFstar(:,:,1:d);
Linf  = Linf(:,:,1:d);
Lstar = Lstar(:,:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
notsteady = 1;
while notsteady && t<smpl
    t = t+1;
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));                   % make sure P is symmetric
    di = data_index{t};
    if isempty(di)
        atilde(:,t)     = a(:,t);
        L(:,:,t)        = T;
        P(:,:,t+1)      = T*P(:,:,t)*T' + QQ;                               %p. 111, DK(2012)
    else
        ZZ = Z(di,:);
        v(di,t)      = Y(di,t) - ZZ*a(:,t);
        F = ZZ*P(:,:,t)*ZZ' + H(di,di);
        sig=sqrt(diag(F));

        if any(diag(F)<kalman_tol) || rcond(F./(sig*sig')) < kalman_tol
            alphahat = Inf;
            return
        end
        iF(di,di,t)   = inv(F./(sig*sig'))./(sig*sig');
        PZI         = P(:,:,t)*ZZ'*iF(di,di,t);
        atilde(:,t) = a(:,t) + PZI*v(di,t);
        K(:,di,t)    = T*PZI;
        L(:,:,t)    = T-K(:,di,t)*ZZ;
        P(:,:,t+1)  = T*P(:,:,t)*L(:,:,t)' + QQ;
    end
    a(:,t+1)    = T*atilde(:,t);
    Pf          = P(:,:,t);
    aK(1,:,t+1) = a(:,t+1);
    for jnk=1:nk
        Pf = T*Pf*T' + QQ;
        PK(jnk,:,:,t+jnk) = Pf;
        if jnk>1
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        end
    end
    %    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<kalman_tol);
end
% $$$ if t<smpl
% $$$     PZI_s = PZI;
% $$$     K_s = K(:,:,t);
% $$$     iF_s = iF(:,:,t);
% $$$     P_s = P(:,:,t+1);
% $$$     P  = cat(3,P(:,:,1:t),repmat(P_s,[1 1 smpl-t]));
% $$$     iF = cat(3,iF(:,:,1:t),repmat(iF_s,[1 1 smpl-t]));
% $$$     L  = cat(3,L(:,:,1:t),repmat(T-K_s*Z,[1 1 smpl-t]));
% $$$     K  = cat(3,K(:,:,1:t),repmat(T*P_s*Z'*iF_s,[1 1 smpl-t]));
% $$$ end
% $$$ while t<smpl
% $$$     t=t+1;
% $$$     v(:,t) = Y(:,t) - Z*a(:,t);
% $$$     atilde(:,t) = a(:,t) + PZI*v(:,t);
% $$$     a(:,t+1) = T*atilde(:,t);
% $$$     Pf          = P(:,:,t);
% $$$     for jnk=1:nk,
% $$$   Pf = T*Pf*T' + QQ;
% $$$         aK(jnk,:,t+jnk) = T^jnk*atilde(:,t);
% $$$   PK(jnk,:,:,t+jnk) = Pf;
% $$$     end
% $$$ end
%% backward pass; r_T and N_T, stored in entry (smpl+1) were initialized at 0
t = smpl+1;
while t>d+1
    t = t-1;
    di = data_index{t};
    if isempty(di)
        % in this case, L is simply T due to Z=0, so that DK (2012), eq. 4.93 obtains
        r(:,t) = L(:,:,t)'*r(:,t+1);                                        %compute r_{t-1}, DK (2012), eq. 4.38 with Z=0
        if state_uncertainty_flag
            N(:,:,t)=L(:,:,t)'*N(:,:,t+1)*L(:,:,t);                         %compute N_{t-1}, DK (2012), eq. 4.42 with Z=0
        end
    else
        ZZ = Z(di,:);
        r(:,t) = ZZ'*iF(di,di,t)*v(di,t) + L(:,:,t)'*r(:,t+1);              %compute r_{t-1}, DK (2012), eq. 4.38
        if state_uncertainty_flag
            N(:,:,t)=ZZ'*iF(di,di,t)*ZZ+L(:,:,t)'*N(:,:,t+1)*L(:,:,t);      %compute N_{t-1}, DK (2012), eq. 4.42
        end
    end
    alphahat(:,t)       = a(:,t) + P(:,:,t)*r(:,t);                         %DK (2012), eq. 4.35
    etahat(:,t) = QRt*r(:,t);                                               %DK (2012), eq. 4.63
    if state_uncertainty_flag
        V(:,:,t)    = P(:,:,t)-P(:,:,t)*N(:,:,t)*P(:,:,t);                      %DK (2012), eq. 4.43
    end
end
if d %diffuse periods
     % initialize r_d^(0) and r_d^(1) as below DK (2012), eq. 5.23
    r0 = zeros(mm,d+1);
    r0(:,d+1) = r(:,d+1);   %set r0_{d}, i.e. shifted by one period
    r1 = zeros(mm,d+1);     %set r1_{d}, i.e. shifted by one period
    if state_uncertainty_flag
        %N_0 at (d+1) is N(d+1), so we can use N for continuing and storing N_0-recursion
        N_1=zeros(mm,mm,d+1);   %set N_1_{d}=0, i.e. shifted by one period, below  DK (2012), eq. 5.26
        N_2=zeros(mm,mm,d+1);   %set N_2_{d}=0, i.e. shifted by one period, below  DK (2012), eq. 5.26
    end
    for t = d:-1:1
        di = data_index{t};
        if isempty(di)
            r1(:,t) = Linf(:,:,t)'*r1(:,t+1);
        else
            if ~Finf_singular(1,t)
                r0(:,t) = Linf(:,:,t)'*r0(:,t+1);                                   % DK (2012), eq. 5.21 where L^(0) is named Linf
                r1(:,t) = Z(di,:)'*(iFinf(di,di,t)*v(di,t)-Kstar(:,di,t)'*T'*r0(:,t+1)) ...
                          + Linf(:,:,t)'*r1(:,t+1);                                       % DK (2012), eq. 5.21, noting that i) F^(1)=(F^Inf)^(-1)(see 5.10), ii) where L^(0) is named Linf, and iii) Kstar=T^{-1}*K^(1)
                if state_uncertainty_flag
                    L_1=(-T*Kstar(:,di,t)*Z(di,:));                                     % noting that Kstar=T^{-1}*K^(1)
                    N(:,:,t)=Linf(:,:,t)'*N(:,:,t+1)*Linf(:,:,t);                       % DK (2012), eq. 5.19, noting that L^(0) is named Linf
                    N_1(:,:,t)=Z(di,:)'*iFinf(di,di,t)*Z(di,:)+Linf(:,:,t)'*N_1(:,:,t+1)*Linf(:,:,t)...
                        +L_1'*N(:,:,t+1)*Linf(:,:,t);                                   % DK (2012), eq. 5.29; note that, compared to DK (2003) this drops the term (L_1'*N(:,:,t+1)*Linf(:,:,t))' in the recursion due to it entering premultiplied by Pinf when computing V, and Pinf*Linf'*N=0
                    N_2(:,:,t)=Z(di,:)'*(-iFinf(di,di,t)*Fstar(di,di,t)*iFinf(di,di,t))*Z(di,:) ...
                        + Linf(:,:,t)'*N_2(:,:,t+1)*Linf(:,:,t)...
                        + Linf(:,:,t)'*N_1(:,:,t+1)*L_1...
                        + L_1'*N_1(:,:,t+1)'*Linf(:,:,t)...
                        + L_1'*N(:,:,t+1)*L_1;                            % DK (2012), eq. 5.29
                end
            else
                r0(:,t) = Z(di,:)'*iFstar(di,di,t)*v(di,t)-Lstar(:,di,t)'*r0(:,t+1); % DK (2003), eq. (14)
                r1(:,t) = T'*r1(:,t+1);                                             % DK (2003), eq. (14)
                if state_uncertainty_flag
                    N(:,:,t)=Z(di,:)'*iFstar(di,di,t)*Z(di,:)...
                             +Lstar(:,:,t)'*N(:,:,t+1)*Lstar(:,:,t);                     % DK (2003), eq. (14)
                    N_1(:,:,t)=T'*N_1(:,:,t+1)*Lstar(:,:,t);                            % DK (2003), eq. (14)
                    N_2(:,:,t)=T'*N_2(:,:,t+1)*T';                                      % DK (2003), eq. (14)
                end
            end
        end
        alphahat(:,t)   = a(:,t) + Pstar(:,:,t)*r0(:,t) + Pinf(:,:,t)*r1(:,t);      % DK (2012), eq. 5.23
        etahat(:,t)     = QRt*r0(:,t);                                              % DK (2012), p. 135
        if state_uncertainty_flag
            V(:,:,t)=Pstar(:,:,t)-Pstar(:,:,t)*N(:,:,t)*Pstar(:,:,t)...
                     -(Pinf(:,:,t)*N_1(:,:,t)*Pstar(:,:,t))'...
                     - Pinf(:,:,t)*N_1(:,:,t)*Pstar(:,:,t)...
                     - Pinf(:,:,t)*N_2(:,:,t)*Pinf(:,:,t);                                   % DK (2012), eq. 5.30
        end
    end
end

if decomp_flag
    decomp = zeros(nk,mm,rr,smpl+nk);
    ZRQinv = inv(Z*QQ*Z');
    for t = max(d,1):smpl
        di = data_index{t};
        % calculate eta_tm1t
        eta_tm1t = QRt*Z(di,:)'*iF(di,di,t)*v(di,t);
        AAA = P(:,:,t)*Z(di,:)'*ZRQinv(di,di)*bsxfun(@times,Z(di,:)*R,eta_tm1t');
        % calculate decomposition
        decomp(1,:,:,t+1) = AAA;
        for h = 2:nk
            AAA = T*AAA;
            decomp(h,:,:,t+h) = AAA;
        end
    end
end

epsilonhat = Y-Z*alphahat;
