function [alphahat,epsilonhat,etahat,a,P,aK,PK,decomp,V] = missing_DiffuseKalmanSmootherH3_Z(a_initial,T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,diffuse_kalman_tol,decomp_flag,state_uncertainty_flag)
% function [alphahat,epsilonhat,etahat,a1,P,aK,PK,d,decomp] = missing_DiffuseKalmanSmootherH3_Z(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,decomp_flag,state_uncertainty_flag)
% Computes the diffuse kalman smoother in the case of a singular var-cov matrix.
% Univariate treatment of multivariate time series.
%
% INPUTS
%    a_initial:mm*1 vector of initial states
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
%    state_uncertainty_flag     if true, compute uncertainty about smoothed
%                               state estimate
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
%    V:        3D array of state uncertainty matrices
%
% Notes:
%   Outputs are stored in decision-rule order, i.e. to get variables in order of declaration
%   as in M_.endo_names, ones needs code along the lines of:
%   variables_declaration_order(dr.order_var,:) = alphahat
%
% Algorithm:
% 
%   Uses the univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 6.4 + 7.2.5
%   and 
%   Koopman/Durbin (2000): "Fast Filtering and Smoothing for Multivariatze State Space
%   Models", in Journal of Time Series Analysis, vol. 21(3), pp. 281-296.
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
a(:,1)          = a_initial;
a1(:,1)         = a_initial;
aK              = zeros(nk,mm,smpl+nk);

Fstar           = zeros(pp,smpl);
Finf            = zeros(pp,smpl);
Fi              = zeros(pp,smpl);
Ki              = zeros(mm,pp,smpl);
Kstar           = zeros(mm,pp,smpl);
Kinf            = zeros(spstar(1),pp,smpl);
P               = zeros(mm,mm,smpl+1);
P1              = P;
PK              = zeros(nk,mm,mm,smpl+nk);
Pstar           = zeros(spstar(1),spstar(2),smpl); 
Pstar(:,:,1)    = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl); 
Pinf(:,:,1)     = Pinf1;
Pstar1          = Pstar;
Pinf1           = Pinf;
rr              = size(Q,1); % number of structural shocks
QQ              = R*Q*transpose(R);
QRt             = Q*transpose(R);
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
epsilonhat      = zeros(rr,smpl);
r               = zeros(mm,smpl);
if state_uncertainty_flag
    V               = zeros(mm,mm,smpl);
    N               = zeros(mm,mm,smpl);
else
    V=[];
end

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
        v(i,t)      = Y(i,t)-Zi*a(:,t);                                     % nu_{t,i} in 6.13 in DK (2012) 
        Fstar(i,t)  = Zi*Pstar(:,:,t)*Zi' +H(i);                            % F_{*,t} in 5.7 in DK (2012), relies on H being diagonal
        Finf(i,t)   = Zi*Pinf(:,:,t)*Zi';                                   % F_{\infty,t} in 5.7 in DK (2012)
        Kstar(:,i,t) = Pstar(:,:,t)*Zi';                                    % KD (2000), eq. (15)
        if Finf(i,t) > diffuse_kalman_tol && newRank                        % F_{\infty,t,i} = 0, use upper part of bracket on p. 175 DK (2012) for w_{t,i}
            icc=icc+1;
            Kinf(:,i,t)       = Pinf(:,:,t)*Zi';                            % KD (2000), eq. (15)                         
            Kinf_Finf         = Kinf(:,i,t)/Finf(i,t);
            a(:,t)            = a(:,t) + Kinf_Finf*v(i,t);                  % KD (2000), eq. (16)
            Pstar(:,:,t)      = Pstar(:,:,t) + ...
                Kinf(:,i,t)*Kinf_Finf'*(Fstar(i,t)/Finf(i,t)) - ...
                Kstar(:,i,t)*Kinf_Finf' - ...
                Kinf_Finf*Kstar(:,i,t)';                                    % KD (2000), eq. (16)
            Pinf(:,:,t)       = Pinf(:,:,t) - Kinf(:,i,t)*Kinf(:,i,t)'/Finf(i,t); % KD (2000), eq. (16)
        elseif Fstar(i,t) > kalman_tol 
            a(:,t)            = a(:,t) + Kstar(:,i,t)*v(i,t)/Fstar(i,t);    % KD (2000), eq. (17)
            Pstar(:,:,t)      = Pstar(:,:,t) - Kstar(:,i,t)*Kstar(:,i,t)'/Fstar(i,t);   % KD (2000), eq. (17)
            % Pinf is passed through unaltered, see eq. (17) of
            % Koopman/Durbin (2000)
        else
            % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
            % p. 157, DK (2012)        
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
    if newRank
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
        v(i,t)  = Y(i,t) - Zi*a(:,t);                                       % nu_{t,i} in 6.13 in DK (2012)
        Fi(i,t) = Zi*P(:,:,t)*Zi' + H(i);                                   % F_{t,i} in 6.13 in DK (2012), relies on H being diagonal 
        Ki(:,i,t) = P(:,:,t)*Zi';                                           % K_{t,i}*F_(i,t) in 6.13 in DK (2012)
        if Fi(i,t) > kalman_tol
            a(:,t) = a(:,t) + Ki(:,i,t)*v(i,t)/Fi(i,t);                     %filtering according to (6.13) in DK (2012)
            P(:,:,t) = P(:,:,t) - Ki(:,i,t)*Ki(:,i,t)'/Fi(i,t);             %filtering according to (6.13) in DK (2012)
        else
            % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
            % p. 157, DK (2012)
        end
    end
    a1(:,t+1) = T*a(:,t);                                                   %transition according to (6.14) in DK (2012) 
    Pf          = P(:,:,t);
    aK(1,:,t+1) = a1(:,t+1); 
    for jnk=1:nk
        Pf = T*Pf*T' + QQ;
        PK(jnk,:,:,t+jnk) = Pf;
        if jnk>1
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));    
        end
    end
    P(:,:,t+1) = T*P(:,:,t)*T' + QQ;                                        %transition according to (6.14) in DK (2012) 
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

%% do backward pass
ri=zeros(mm,1);
if state_uncertainty_flag
    Ni=zeros(mm,mm);
end
t = smpl+1;
while t > d+1
    t = t-1;
    di = flipud(data_index{t})';
    for i = di
        if Fi(i,t) > kalman_tol
            Li = eye(mm)-Ki(:,i,t)*Z(i,:)/Fi(i,t);
            ri = Z(i,:)'/Fi(i,t)*v(i,t)+Li'*ri;                             % DK (2012), 6.15, equation for r_{t,i-1}
            if state_uncertainty_flag
                Ni = Z(i,:)'/Fi(i,t)*Z(i,:)+Li'*Ni*Li;                      % KD (2000), eq. (23)
            end
        end
    end
    r(:,t) = ri;                                                            % DK (2012), below 6.15, r_{t-1}=r_{t,0} 
    alphahat(:,t) = a1(:,t) + P1(:,:,t)*r(:,t);
    etahat(:,t) = QRt*r(:,t);
    ri = T'*ri;                                                             % KD (2003), eq. (23), equation for r_{t-1,p_{t-1}} 
    if state_uncertainty_flag
        N(:,:,t) = Ni;                                                          % DK (2012), below 6.15, N_{t-1}=N_{t,0} 
        V(:,:,t) = P1(:,:,t)-P1(:,:,t)*N(:,:,t)*P1(:,:,t);                      % KD (2000), eq. (7) with N_{t-1} stored in N(:,:,t)
        Ni = T'*Ni*T;                                                           % KD (2000), eq. (23), equation for N_{t-1,p_{t-1}} 
    end
end
if d
    r0 = zeros(mm,d); 
    r0(:,d) = ri;
    r1 = zeros(mm,d);
    if state_uncertainty_flag
        %N_0 at (d+1) is N(d+1), so we can use N for continuing and storing N_0-recursion
        N_0=zeros(mm,mm,d);   %set N_1_{d}=0, below  KD (2000), eq. (24)
        N_0(:,:,d) = Ni;
        N_1=zeros(mm,mm,d);   %set N_1_{d}=0, below  KD (2000), eq. (24)
        N_2=zeros(mm,mm,d);   %set N_2_{d}=0, below  KD (2000), eq. (24)
    end
    for t = d:-1:1
        di = flipud(data_index{t})';
        for i = di
            if Finf(i,t) > diffuse_kalman_tol
                % recursions need to be from highest to lowest term in order to not 
                % overwrite lower terms still needed in this step
                Linf    = eye(mm) - Kinf(:,i,t)*Z(i,:)/Finf(i,t);
                L0      = (Kinf(:,i,t)*(Fstar(i,t)/Finf(i,t))-Kstar(:,i,t))*Z(i,:)/Finf(i,t);            
                r1(:,t) = Z(i,:)'*v(i,t)/Finf(i,t) + ...
                          L0'*r0(:,t) + ...
                          Linf'*r1(:,t);   % KD (2000), eq. (25) for r_1
                r0(:,t) = Linf'*r0(:,t);   % KD (2000), eq. (25) for r_0
                if state_uncertainty_flag
                    N_2(:,:,t)=Z(i,:)'/Finf(i,t)^2*Z(i,:)*Fstar(i,t) ...
                        + Linf'*N_2(:,:,t)*Linf...
                        + Linf'*N_1(:,:,t)*L0...
                        + L0'*N_1(:,:,t)'*Linf...
                        + L0'*N_0(:,:,t)*L0;                                    % DK (2012), eq. 5.29
                    N_1(:,:,t)=Z(i,:)'/Finf(i,t)*Z(i,:)+Linf'*N_1(:,:,t)*Linf...
                        +L0'*N_0(:,:,t)*Linf;                                   % DK (2012), eq. 5.29; note that, compared to DK (2003) this drops the term (L_1'*N(:,:,t+1)*Linf(:,:,t))' in the recursion due to it entering premultiplied by Pinf when computing V, and Pinf*Linf'*N=0
                    N_0(:,:,t)=Linf'*N_0(:,:,t)*Linf;                           % DK (2012), eq. 5.19, noting that L^(0) is named Linf
                end
            elseif Fstar(i,t) > kalman_tol % step needed whe Finf == 0
                L_i=eye(mm) - Kstar(:,i,t)*Z(i,:)/Fstar(i,t);
                r0(:,t) = Z(i,:)'/Fstar(i,t)*v(i,t)+L_i'*r0(:,t);           % propagate r0 and keep r1 fixed
                if state_uncertainty_flag
                    N_0(:,:,t)=Z(i,:)'/Fstar(i,t)*Z(i,:)+L_i'*N_0(:,:,t)*L_i;   % propagate N_0 and keep N_1 and N_2 fixed
                end
            end
        end
        alphahat(:,t) = a1(:,t) + Pstar1(:,:,t)*r0(:,t) + Pinf1(:,:,t)*r1(:,t); % KD (2000), eq. (26)
        r(:,t)        = r0(:,t);
        etahat(:,t)   = QRt*r(:,t);                                         % KD (2000), eq. (27)
        if state_uncertainty_flag
            V(:,:,t)=Pstar(:,:,t)-Pstar(:,:,t)*N_0(:,:,t)*Pstar(:,:,t)...
                -(Pinf(:,:,t)*N_1(:,:,t)*Pstar(:,:,t))'...
                - Pinf(:,:,t)*N_1(:,:,t)*Pstar(:,:,t)...
                - Pinf(:,:,t)*N_2(:,:,t)*Pinf(:,:,t);                       % DK (2012), eq. 5.30
        end
        if t > 1
            r0(:,t-1) = T'*r0(:,t);                                         % KD (2000), below eq. (25) r_{t-1,p_{t-1}}=T'*r_{t,0}
            r1(:,t-1) = T'*r1(:,t);                                         % KD (2000), below eq. (25) r_{t-1,p_{t-1}}=T'*r_{t,0}
            if state_uncertainty_flag
                N_0(:,:,t-1)= T'*N_0(:,t)*T;                                % KD (2000), below eq. (25) N_{t-1,p_{t-1}}=T'*N_{t,0}*T
                N_1(:,:,t-1)= T'*N_1(:,t)*T;                                % KD (2000), below eq. (25) N^1_{t-1,p_{t-1}}=T'*N^1_{t,0}*T
                N_2(:,:,t-1)= T'*N_2(:,t)*T;                                % KD (2000), below eq. (25) N^2_{t-1,p_{t-1}}=T'*N^2_{t,0}*T
            end
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
