function [dLIK,dlik,a,Pstar] = missing_observations_kalman_filter_d(data_index,number_of_observations,no_more_missing_observations, ...
                                                      Y, start, last, ...
                                                      a, Pinf, Pstar, ...
                                                      kalman_tol, riccati_tol, presample, ...
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
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98). 

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


% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
dlik = zeros(smpl,1);      % Initialization of the vector gathering the densities.
dLIK = Inf;                % Default value of the log likelihood.
oldK = Inf;
crit1=1.e-6;

if isequal(H,0)
    H = zeros(pp,pp);
end
s = 0;

while rank(Pinf,crit1) && (t<=last)
    s = t-start+1;
    d_index = data_index{t};
    if isempty(d_index)
        a = T*a;
        Pstar = T*Pstar*transpose(T)+QQ;
        Pinf  = T*Pinf*transpose(T);
    else
        ZZ = Z(d_index,:);
        v  = Y(d_index,t)-ZZ*a;
        Finf  = ZZ*Pinf*ZZ';
        if rcond(Finf) < kalman_tol
            if ~all(abs(Finf(:)) < kalman_tol)
                % The univariate diffuse kalman filter should be used.
                return
            else
                Fstar = ZZ*Pstar*ZZ' + H(d_index,d_index);
                if rcond(Fstar) < kalman_tol
                    if ~all(abs(Fstar(:))<kalman_tol)
                        % The univariate diffuse kalman filter should be used.
                        return
                    else
                        a = T*a;
                        Pstar = T*Pstar*transpose(T)+QQ;
                        Pinf  = T*Pinf*transpose(T);
                    end
                else
                    iFstar = inv(Fstar);
                    dFstar = det(Fstar);
                    Kstar  = Pstar*ZZ'*iFstar;
                    dlik(s) = log(dFstar) + v'*iFstar*v + length(d_index)*log(2*pi);
                    Pinf   = T*Pinf*transpose(T);
                    Pstar  = T*(Pstar-Pstar*ZZ'*Kstar')*T'+QQ;
                    a      = T*(a+Kstar*v);
                end
            end
        else
            dlik(s) = log(det(Finf));
            iFinf  = inv(Finf);
            Kinf   = Pinf*ZZ'*iFinf;
            Fstar  = ZZ*Pstar*ZZ' + H(d_index,d_index);
            Kstar  = (Pstar*ZZ'-Kinf*Fstar)*iFinf;
            Pstar  = T*(Pstar-Pstar*ZZ'*Kinf'-Pinf*ZZ'*Kstar')*T'+QQ;
            Pinf   = T*(Pinf-Pinf*ZZ'*Kinf')*T';
            a      = T*(a+Kinf*v);
        end
    end
    t  = t+1;
end

if t==(last+1)
    warning(['There isn''t enough information to estimate the initial conditions of the nonstationary variables']);                   
    dLIK = NaN;
    return
end

dlik = .5*dlik(1:s);

dLIK = sum(dlik(1+presample:end));
