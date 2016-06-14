function [LIK, lik,a,P] = univariate_kalman_filter(data_index,number_of_observations,no_more_missing_observations,Y,start,last,a,P,kalman_tol,riccati_tol,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods,analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P)
% Computes the likelihood of a stationnary state space model (univariate approach).

%@info:
%! @deftypefn {Function File} {[@var{LIK},@var{likk},@var{a},@var{P} ] =} univariate_kalman_filter (@var{data_index}, @var{number_of_observations},@var{no_more_missing_observations}, @var{Y}, @var{start}, @var{last}, @var{a}, @var{P}, @var{kalman_tol}, @var{riccati_tol},@var{presample},@var{T},@var{Q},@var{R},@var{H},@var{Z},@var{mm},@var{pp},@var{rr},@var{Zflag},@var{diffuse_periods})
%! @anchor{univariate_kalman_filter}
%! @sp 1
%! Computes the likelihood of a stationary state space model, given initial condition for the states (mean and variance).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item data_index
%! Matlab's cell, 1*T cell of column vectors of indices (in the vector of observed variables).
%! @item number_of_observations
%! Integer scalar, effective number of observations.
%! @item no_more_missing_observations
%! Integer scalar, date after which there is no more missing observation (it is then possible to switch to the steady state kalman filter).
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
%! Vector (@var{pp}) of doubles, diagonal of covariance matrix of the measurement errors (corelation among measurement errors is handled by a model transformation).
%! @item Z
%! Matrix (@var{pp}*@var{mm}) of doubles or vector of integers, matrix relating the states to the observed variables or vector of indices (depending on the value of @var{Zflag}).
%! @item mm
%! Integer scalar, number of state variables.
%! @item pp
%! Integer scalar, number of observed variables.
%! @item rr
%! Integer scalar, number of structural innovations.
%! @item Zflag
%! Integer scalar, equal to 0 if Z is a vector of indices targeting the observed variables in the state vector, equal to 1 if Z is a @var{pp}*@var{mm} matrix.
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
%! @ref{dsge_likelihood}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{univariate_kalman_filter_ss}
%! @end deftypefn
%@eod:
% 
% Algorithm:
% 
%   Uses the univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 6.4 + 7.2.5


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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

if nargin<20 || isempty(Zflag)% Set default value for Zflag ==> Z is a vector of indices.
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<21
    diffuse_periods = 0;
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
lik  = zeros(smpl,pp);     % Initialization of the matrix gathering the densities at each time and each observable
LIK  = Inf;                % Default value of the log likelihood.
oldP = Inf;
l2pi = log(2*pi);
notsteady = 1;

oldK = Inf;
K = NaN(mm,pp);
asy_hess=0;

if  analytic_derivation == 0,
    DLIK=[];
    Hess=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zeros(k,1);                             % Initialization of the score.
    Da    = zeros(mm,k);                            % Derivative State vector.
    dlik  = zeros(smpl,k);
    
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
end

while notsteady && t<=last %loop over t
    s = t-start+1;
    d_index = data_index{t};
    if Zflag
        z = Z(d_index,:);
    else
        z = Z(d_index);
    end
    oldP = P(:);
    for i=1:rows(z) %loop over i
        if Zflag
            prediction_error = Y(d_index(i),t) - z(i,:)*a;  % nu_{t,i} in 6.13 in DK (2012) 
            PZ = P*z(i,:)';                                 % Z_{t,i}*P_{t,i}*Z_{t,i}'
            Fi = z(i,:)*PZ + H(d_index(i));                 % F_{t,i} in 6.13 in DK (2012), relies on H being diagonal 
        else
            prediction_error = Y(d_index(i),t) - a(z(i));   % nu_{t,i} in 6.13 in DK (2012) 
            PZ = P(:,z(i));                                 % Z_{t,i}*P_{t,i}*Z_{t,i}'
            Fi = PZ(z(i)) + H(d_index(i));                  % F_{t,i} in 6.13 in DK (2012), relies on H being diagonal 
        end
        if Fi>kalman_tol
            Ki =  PZ/Fi; %K_{t,i} in 6.13 in DK (2012)
            if t>=no_more_missing_observations
                K(:,i) = Ki;
            end
            lik(s,i) = log(Fi) + (prediction_error*prediction_error)/Fi + l2pi; %Top equation p. 175 in DK (2012)
            if analytic_derivation,
                if analytic_derivation==2,
                    [Da,DP,DLIKt,D2a,D2P, Hesst] = univariate_computeDLIK(k,i,z(i,:),Zflag,prediction_error,Ki,PZ,Fi,Da,DYss,DP,DH(d_index(i),:),notsteady,D2a,D2Yss,D2P);
                else
                    [Da,DP,DLIKt,Hesst] = univariate_computeDLIK(k,i,z(i,:),Zflag,prediction_error,Ki,PZ,Fi,Da,DYss,DP,DH(d_index(i),:),notsteady);
                end
                if t>presample
                    DLIK = DLIK + DLIKt;
                    if analytic_derivation==2 || asy_hess,
                        Hess = Hess + Hesst;
                    end
                end
                dlik(s,:)=dlik(s,:)+DLIKt';
            end
            a = a + Ki*prediction_error; %filtering according to (6.13) in DK (2012)
            P = P - PZ*Ki';              %filtering according to (6.13) in DK (2012)
        else
            % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
            % p. 157, DK (2012)
        end
    end
    if analytic_derivation,        
        if analytic_derivation==2,
            [Da,DP,D2a,D2P] = univariate_computeDstate(k,a,P,T,Da,DP,DT,DOm,notsteady,D2a,D2P,D2T,D2Om);
        else
            [Da,DP] = univariate_computeDstate(k,a,P,T,Da,DP,DT,DOm,notsteady);
        end
    end
    a = T*a;            %transition according to (6.14) in DK (2012) 
    P = T*P*T' + QQ;    %transition according to (6.14) in DK (2012)
    if t>=no_more_missing_observations
        notsteady = max(abs(K(:)-oldK))>riccati_tol;
        oldK = K(:);
    end
    t = t+1;
end

% Divide by two.
lik(1:s,:) = .5*lik(1:s,:);
if analytic_derivation,
    DLIK = DLIK/2;
    dlik = dlik/2;
    if analytic_derivation==2 || asy_hess,
%         Hess = (Hess + Hess')/2;
        Hess = -Hess/2;
    end
end

% Call steady state univariate kalman filter if needed.
if t <= last
    if analytic_derivation,
        if analytic_derivation==2,
            [tmp, tmp2] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag, ...
                analytic_derivation,Da,DT,DYss,DP,DH,D2a,D2T,D2Yss,D2P);
        else
            [tmp, tmp2] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag, ...
                analytic_derivation,Da,DT,DYss,DP,DH,asy_hess);
        end
        lik(s+1:end,:)=tmp2{1};
        dlik(s+1:end,:)=tmp2{2};
        DLIK = DLIK + tmp{2};
        if analytic_derivation==2 || asy_hess,
            Hess = Hess + tmp{3};
        end
    else
        [tmp, lik(s+1:end,:)] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag);
    end
end

% Compute minus the log-likelihood.
if presample > diffuse_periods
    LIK = sum(sum(lik(1+presample-diffuse_periods:end,:)));
else
    LIK = sum(sum(lik));
end

if analytic_derivation,
    if analytic_derivation==2 || asy_hess,
        LIK={LIK, DLIK, Hess};
    else
        LIK={LIK, DLIK};
    end
    lik={lik, dlik};
end
