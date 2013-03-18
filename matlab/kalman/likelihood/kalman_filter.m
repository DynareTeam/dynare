function [LIK, LIKK, a, P] = kalman_filter(Y,start,last,a,P,kalman_tol,riccati_tol,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods,analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P)
% Computes the likelihood of a stationary state space model.

%@info:
%! @deftypefn {Function File} {[@var{LIK},@var{likk},@var{a},@var{P} ] =} DsgeLikelihood (@var{Y}, @var{start}, @var{last}, @var{a}, @var{P}, @var{kalman_tol}, @var{riccati_tol},@var{presample},@var{T},@var{Q},@var{R},@var{H},@var{Z},@var{mm},@var{pp},@var{rr},@var{Zflag},@var{diffuse_periods})
%! @anchor{kalman_filter}
%! @sp 1
%! Computes the likelihood of a stationary state space model, given initial condition for the states (mean and variance).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item Y
%! Matrix (@var{pp}*T) of doubles, data.
%! @item start
%! Integer scalar, first period.
%! @item last
%! Integer scalar, last period (@var{last}-@var{first} has to be inferior to T).
%! @item a
%! Vector (@var{mm}*1) of doubles, initial mean of the state vector.
%! @item P
%! Matrix (@var{mm}*@var{mm}) of doubles, initial covariance matrix of the state vector.
%! @item kalman_tol
%! Double scalar, tolerance parameter (rcond, inversibility of the covariance matrix of the prediction errors).
%! @item riccati_tol
%! Double scalar, tolerance parameter (iteration over the Riccati equation).
%! @item presample
%! Integer scalar, presampling if strictly positive (number of initial iterations to be discarded when evaluating the likelihood).
%! @item T
%! Matrix (@var{mm}*@var{mm}) of doubles, transition matrix of the state equation.
%! @item Q
%! Matrix (@var{rr}*@var{rr}) of doubles, covariance matrix of the structural innovations (noise in the state equation).
%! @item R
%! Matrix (@var{mm}*@var{rr}) of doubles, second matrix of the state equation relating the structural innovations to the state variables.
%! @item H
%! Matrix (@var{pp}*@var{pp}) of doubles, covariance matrix of the measurement errors (if no measurement errors set H as a zero scalar).
%! @item Z
%! Matrix (@var{pp}*@var{mm}) of doubles or vector of integers, matrix relating the states to the observed variables or vector of indices (depending on the value of @var{Zflag}).
%! @item mm
%! Integer scalar, number of state variables.
%! @item pp
%! Integer scalar, number of observed variables.
%! @item rr
%! Integer scalar, number of structural innovations.
%! @item Zflag
%! Integer scalar, equal to 0 if Z is a vector of indices targeting the obseved variables in the state vector, equal to 1 if Z is a @var{pp}*@var{mm} matrix.
%! @item diffuse_periods
%! Integer scalar, number of diffuse filter periods in the initialization step.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item LIK
%! Double scalar, value of (minus) the likelihood.
%! @item likk
%! Column vector of doubles, values of the density of each observation.
%! @item a
%! Vector (@var{mm}*1) of doubles, mean of the state vector at the end of the (sub)sample.
%! @item P
%! Matrix (@var{mm}*@var{mm}) of doubles, covariance of the state vector at the end of the (sub)sample.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{DsgeLikelihood}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{kalman_filter_ss}
%! @end deftypefn
%@eod:

% Copyright (C) 2004-2012 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

% Set defaults.
if nargin<17
    Zflag = 0;
end

if nargin<18
    diffuse_periods = 0;
end

if nargin<19
    analytic_derivation = 0;
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
likk = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;
F_singular  = 1;
asy_hess=0;

if  analytic_derivation == 0,
    DLIK=[];
    Hess=[];
    LIKK=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zeros(k,1);                             % Initialization of the score.
    Da    = zeros(mm,k);                            % Derivative State vector.
    dlikk = zeros(smpl,k);
    
    if Zflag==0,
        C = zeros(pp,mm);
        for ii=1:pp; C(ii,Z(ii))=1;end         % SELECTION MATRIX IN MEASUREMENT EQ. (FOR WHEN IT IS NOT CONSTANT)
    else
        C=Z;
    end
    dC = zeros(pp,mm,k);   % either selection matrix or schur have zero derivatives
    if analytic_derivation==2,
        Hess  = zeros(k,k);                             % Initialization of the Hessian
        D2a    = zeros(mm,k,k);                             % State vector.
        d2C = zeros(pp,mm,k,k);
    else
        asy_hess=D2T;
        Hess=[];
        D2a=[];
        D2T=[];
        D2Yss=[];
    end
    if asy_hess,
        Hess  = zeros(k,k);                             % Initialization of the Hessian
    end
    LIK={inf,DLIK,Hess};
    LIKK={likk,dlikk};
end

while notsteady && t<=last
    s = t-start+1;
    if Zflag
        v  = Y(:,t)-Z*a;
        F  = Z*P*Z' + H;
    else
        v  = Y(:,t)-a(Z);
        F  = P(Z,Z) + H;
    end
    if rcond(F) < kalman_tol
        if ~all(abs(F(:))<kalman_tol)
            return
        else
            a = T*a;
            P = T*P*transpose(T)+QQ;
        end
    else
        F_singular = 0;
        dF      = det(F);
        iF      = inv(F);
        likk(s) = log(dF)+transpose(v)*iF*v;
        if Zflag
            K = P*Z'*iF;
            Ptmp = T*(P-K*Z*P)*transpose(T)+QQ;
        else
            K = P(:,Z)*iF;
            Ptmp = T*(P-K*P(Z,:))*transpose(T)+QQ;
        end
        tmp = (a+K*v);
        if analytic_derivation,
            if analytic_derivation==2,
                [Da,DP,DLIKt,D2a,D2P, Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2P);
            else
                [Da,DP,DLIKt,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady);
            end
            if t>presample
                DLIK = DLIK + DLIKt;
                if analytic_derivation==2 || asy_hess,
                    Hess = Hess + Hesst;
                end
            end
            dlikk(s,:)=DLIKt;
        end
        a = T*tmp;
        P=Ptmp;
        notsteady = max(abs(K(:)-oldK))>riccati_tol;
        oldK = K(:);
    end
    t = t+1;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Add observation's densities constants and divide by two.
likk(1:s) = .5*(likk(1:s) + pp*log(2*pi));
if analytic_derivation,
    DLIK = DLIK/2;
    dlikk = dlikk/2;
    if analytic_derivation==2 || asy_hess,
        if asy_hess==0,
        Hess = Hess + tril(Hess,-1)';
        end
        Hess = -Hess/2;
    end
end

% Call steady state Kalman filter if needed.
if t <= last
    if analytic_derivation,
        if analytic_derivation==2,
            [tmp, tmp2] = kalman_filter_ss(Y,t,last,a,T,K,iF,dF,Z,pp,Zflag, ...
                analytic_derivation,Da,DT,DYss,D2a,D2T,D2Yss);
        else
            [tmp, tmp2] = kalman_filter_ss(Y,t,last,a,T,K,iF,dF,Z,pp,Zflag, ...
                analytic_derivation,Da,DT,DYss,asy_hess);
        end
        likk(s+1:end)=tmp2{1};
        dlikk(s+1:end,:)=tmp2{2};
        DLIK = DLIK + tmp{2};
        if analytic_derivation==2 || asy_hess,
            Hess = Hess + tmp{3};
        end
    else
        [tmp, likk(s+1:end)] = kalman_filter_ss(Y,t,last,a,T,K,iF,dF,Z,pp,Zflag);
    end
end

% Compute minus the log-likelihood.
if presample>diffuse_periods,
    LIK = sum(likk(1+(presample-diffuse_periods):end));
else
    LIK = sum(likk);
end

if analytic_derivation,
    if analytic_derivation==2 || asy_hess,
        LIK={LIK, DLIK, Hess};
    else
        LIK={LIK, DLIK};
    end
    LIKK={likk, dlikk};
else
    LIKK=likk;
end
