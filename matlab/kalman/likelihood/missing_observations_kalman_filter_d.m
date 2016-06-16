function [dLIK,dlik,a,Pstar] = missing_observations_kalman_filter_d(data_index,number_of_observations,no_more_missing_observations, ...
                                                      Y, start, last, ...
                                                      a, Pinf, Pstar, ...
                                                      kalman_tol, diffuse_kalman_tol, riccati_tol, presample, ...
                                                      T, R, Q, H, Z, mm, pp, rr)
% Computes the diffuse likelihood of a state space model when some observations are missing.
%
% INPUTS 
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%    Y                            [double]      pp*smpl matrix of (detrended) data, where pp is the number of observed variables.
%    start                        [integer]     scalar, first observation.
%    last                         [integer]     scalar, last observation.
%    a                            [double]      mm*1 vector, levels of the state variables.
%    Pinf                         [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    Pstar                        [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    kalman_tol                   [double]      scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]      scalar, tolerance parameter (riccati iteration).
%    presample                    [integer]     scalar, presampling if strictly positive.    
%    T                            [double]      mm*mm matrix, transition matrix in  the state equations.
%    R                            [double]      mm*rr matrix relating the structural innovations to the state vector.
%    Q                            [double]      rr*rr covariance matrix of the structural innovations.
%    H                            [double]      pp*pp covariance matrix of the measurement errors (if H is equal to zero (scalar) there is no measurement error). 
%    Z                            [double]      pp*mm matrix, selection matrix or pp linear independant combinations of the state vector.
%    mm                           [integer]     scalar, number of state variables.
%    pp                           [integer]     scalar, number of observed variables.
%    rr                           [integer]     scalar, number of structural innovations.    
%             
% OUTPUTS 
%    dLIK        [double]    scalar, MINUS loglikelihood
%    dlik        [double]    vector, density of observations in each period.
%    a           [double]    mm*1 vector, estimated level of the states.
%    Pstar       [double]    mm*mm matrix, covariance matrix of the states.        
%        
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003), in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98. 
%   and
%   Durbin/Koopman (2012): "Time Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 5 and 7.2

% 
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


% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
dlik = zeros(smpl,1);      % Initialization of the vector gathering the densities.
dLIK = Inf;                % Default value of the log likelihood.
oldK = Inf;

if isequal(H,0)
    H = zeros(pp,pp);
end
s = 0;

while rank(Pinf,diffuse_kalman_tol) && (t<=last)
    s = t-start+1;
    d_index = data_index{t};
    if isempty(d_index)
        %no observations, propagate forward without updating based on
        %observations
        a = T*a;
        Pstar = T*Pstar*transpose(T)+QQ;
        Pinf  = T*Pinf*transpose(T);
    else
        ZZ = Z(d_index,:);                                                  %span selector matrix
        v  = Y(d_index,t)-ZZ*a;                                             %get prediction error v^(0) in (5.13) DK (2012)         
        Finf  = ZZ*Pinf*ZZ';                                                % (5.7) in DK (2012)
        if rcond(Finf) < diffuse_kalman_tol                                 %F_{\infty,t} = 0 
            if ~all(abs(Finf(:)) < diffuse_kalman_tol)                      %rank-deficient but not rank 0
                % The univariate diffuse kalman filter should be used.
                return
            else                                                            %rank of F_{\infty,t} is 0
                Fstar = ZZ*Pstar*ZZ' + H(d_index,d_index);                  % (5.7) in DK (2012)
                if rcond(Fstar) < kalman_tol                                %F_{*} is singular
                    if ~all(abs(Fstar(:))<kalman_tol)
                        % The univariate diffuse kalman filter should be used.
                        return
                    else %rank 0
                        a = T*a;
                        Pstar = T*Pstar*transpose(T)+QQ;
                        Pinf  = T*Pinf*transpose(T);                        % (5.16) DK (2012)
                    end
                else
                    iFstar = inv(Fstar);
                    dFstar = det(Fstar);
                    Kstar  = Pstar*ZZ'*iFstar;                              %(5.15) of DK (2012) with Kstar=K^(0)*T^{-1}
                    dlik(s) = log(dFstar) + v'*iFstar*v + length(d_index)*log(2*pi);    %set w_t to bottom case in bottom equation page 172, DK (2012)
                    Pinf   = T*Pinf*transpose(T);                           % (5.16) DK (2012)
                    Pstar  = T*(Pstar-Pstar*ZZ'*Kstar')*T'+QQ;              % (5.17) DK (2012)
                    a      = T*(a+Kstar*v);                                 % (5.13) DK (2012)
                end
            end
        else
            dlik(s) = log(det(Finf))+length(d_index)*log(2*pi);             %set w_t to top case in bottom equation page 172, DK (2012)
            iFinf  = inv(Finf);
            Kinf   = Pinf*ZZ'*iFinf;
            %see notes in kalman_filter_d.m for details of computations
            Fstar  = ZZ*Pstar*ZZ' + H(d_index,d_index);                     %(5.7) DK(2012)
            Kstar  = (Pstar*ZZ'-Kinf*Fstar)*iFinf;                          %(5.12) DK(2012)
            Pstar  = T*(Pstar-Pstar*ZZ'*Kinf'-Pinf*ZZ'*Kstar')*T'+QQ;       %(5.14) DK(2012)
            Pinf   = T*(Pinf-Pinf*ZZ'*Kinf')*T';                            %(5.14) DK(2012)
            a      = T*(a+Kinf*v);                                          %(5.13) DK(2012)
        end
    end
    t  = t+1;
end

if t==(last+1)
    warning(['kalman_filter_d: There isn''t enough information to estimate the initial conditions of the nonstationary variables. The diffuse Kalman filter never left the diffuse stage.']);                   
    dLIK = NaN;
    return
end

dlik = .5*dlik(1:s);

dLIK = sum(dlik(1+presample:end));
