function [LIK,likk,a] = univariate_kalman_filter_ss(Y,start,last,a,P,kalman_tol,T,H,Z,pp,Zflag,analytic_derivation,Da,DT,DYss,DP,DH,D2a,D2T,D2Yss,D2P)
% Computes the likelihood of a stationnary state space model (steady state univariate kalman filter).

%@info:
%! @deftypefn {Function File} {[@var{LIK},@var{likk},@var{a} ] =} univariate_kalman_filter_ss (@var{Y}, @var{start}, @var{last}, @var{a}, @var{P}, @var{kalman_tol}, @var{riccati_tol},@var{presample},@var{T},@var{Q},@var{R},@var{H},@var{Z},@var{mm},@var{pp},@var{rr},@var{Zflag},@var{diffuse_periods})
%! @anchor{univariate_kalman_filter_ss}
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
%! Matrix (@var{mm}*@var{mm}) of doubles, steady state covariance matrix of the state vector.
%! @item kalman_tol
%! Double scalar, tolerance parameter (rcond, inversibility of the covariance matrix of the prediction errors).
%! @item T
%! Matrix (@var{mm}*@var{mm}) of doubles, transition matrix of the state equation.
%! @item H
%! Vector (@var{pp}) of doubles, diagonal of covariance matrix of the measurement errors (corelation among measurement errors is handled by a model transformation).
%! Matrix (@var{pp}*@var{pp}) of doubles, covariance matrix of the measurement errors (if no measurement errors set H as a zero scalar).
%! @item Z
%! Matrix (@var{pp}*@var{mm}) of doubles or vector of integers, matrix relating the states to the observed variables or vector of indices (depending on the value of @var{Zflag}).
%! @item pp
%! Integer scalar, number of observed variables.
%! @item Zflag
%! Integer scalar, equal to 0 if Z is a vector of indices targeting the obseved variables in the state vector, equal to 1 if Z is a @var{pp}*@var{mm} matrix.
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
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{univariate_kalman_filter}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @end deftypefn
%@eod:
% 
% Algorithm: See univariate_kalman_filter.m

% Copyright (C) 2011-2016 Dynare Team
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

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
t    = start;              % Initialization of the time index.
likk = zeros(smpl,pp);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
l2pi = log(2*pi);
asy_hess=0;

if nargin<12
    analytic_derivation = 0;
end

if  analytic_derivation == 0,
    DLIK=[];
    Hess=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zeros(k,1);                             % Initialization of the score.
    dlikk = zeros(smpl,k);
    if analytic_derivation==2,
        Hess  = zeros(k,k);                             % Initialization of the Hessian
    else
        asy_hess=D2a;
        if asy_hess,
            Hess  = zeros(k,k);                             % Initialization of the Hessian
        else
            Hess=[];
        end
    end
end

% Steady state kalman filter.
while t<=last
    s  = t-start+1;
    PP = P;
    if analytic_derivation,
        DPP = DP;
        if analytic_derivation==2,
            D2PP = D2P;
        end
    end
    for i=1:pp
        if Zflag
            prediction_error = Y(i,t) - Z(i,:)*a;
            PPZ = PP*Z(i,:)';
            Fi = Z(i,:)*PPZ + H(i);
        else
            prediction_error = Y(i,t) - a(Z(i));
            PPZ = PP(:,Z(i));
            Fi = PPZ(Z(i)) + H(i);
        end
        if Fi>kalman_tol
            Ki = PPZ/Fi;
            a  = a + Ki*prediction_error;
            PP = PP - PPZ*Ki';
            likk(s,i) = log(Fi) + prediction_error*prediction_error/Fi + l2pi;
            if analytic_derivation,
                if analytic_derivation==2,
                    [Da,DPP,DLIKt,D2a,D2PP, Hesst] = univariate_computeDLIK(k,i,Z(i,:),Zflag,prediction_error,Ki,PPZ,Fi,Da,DYss,DPP,DH(i,:),0,D2a,D2Yss,D2PP);
                else
                    [Da,DPP,DLIKt,Hesst] = univariate_computeDLIK(k,i,Z(i,:),Zflag,prediction_error,Ki,PPZ,Fi,Da,DYss,DPP,DH(i,:),0);
                end
                DLIK = DLIK + DLIKt;
                if analytic_derivation==2 || asy_hess,
                    Hess = Hess + Hesst;
                end
                dlikk(s,:)=dlikk(s,:)+DLIKt';
            end
        else
            % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
            % p. 157, DK (2012)
        end
    end
    if analytic_derivation,        
        if analytic_derivation==2,
            [Da,junk,D2a] = univariate_computeDstate(k,a,P,T,Da,DP,DT,[],0,D2a,D2P,D2T);
        else
            Da = univariate_computeDstate(k,a,P,T,Da,DP,DT,[],0);
        end
    end
    a = T*a;
    t = t+1;
end

likk = .5*likk;

LIK = sum(sum(likk));
if analytic_derivation,
    dlikk = dlikk/2;
    DLIK = DLIK/2;
    likk = {likk, dlikk};
end
if analytic_derivation==2 || asy_hess,
%     Hess = (Hess + Hess')/2;
    Hess = -Hess/2;
    LIK={LIK,DLIK,Hess};
elseif analytic_derivation==1,
    LIK={LIK,DLIK};
end