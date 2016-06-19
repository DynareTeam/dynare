function  [LIK, lik, a, P] = missing_observations_kalman_filter(data_index,number_of_observations,no_more_missing_observations,Y,start,last,a,P,kalman_tol,riccati_tol,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods)
% Computes the likelihood of a state space model in the case with missing observations.
%
% INPUTS
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%    Y                            [double]    pp*smpl matrix of data.    
%    start                        [integer]   scalar, index of the first observation.
%    last                         [integer]   scalar, index of the last observation.
%    a                            [double]    pp*1 vector, initial level of the state vector.
%    P                            [double]    pp*pp matrix, covariance matrix of the initial state vector.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).    
%    presample                    [integer]   scalar, presampling if strictly positive.
%    T                            [double]    mm*mm transition matrix of the state equation.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.    
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    H                            [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    Z                            [integer]   pp*1 vector of indices for the observed variables.    
%    mm                           [integer]   scalar, dimension of the state vector.
%    pp                           [integer]   scalar, number of observed variables.
%    rr                           [integer]   scalar, number of structural innovations.    
%    
% OUTPUTS 
%    LIK        [double]    scalar, MINUS loglikelihood
%    lik        [double]    vector, density of observations in each period.
%    a          [double]    mm*1 vector, estimated level of the states.
%    P          [double]    mm*mm matrix, covariance matrix of the states.        
%        
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.   

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

% Set defaults
if nargin<20
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<21
    diffuse_periods = 0;
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end

if isequal(H,0)
    H = zeros(pp,pp);
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
lik  = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;
F_singular  = 1;
s = 0;

while notsteady && t<=last
    s  = t-start+1;
    d_index = data_index{t};
    if isempty(d_index)
        a = T*a;
        P = T*P*transpose(T)+QQ;
    else
        % Compute the prediction error and its variance
        if Zflag
            z = Z(d_index,:);
            v = Y(d_index,t)-z*a;
            F = z*P*z' + H(d_index,d_index);
        else
            z = Z(d_index);
            v = Y(d_index,t) - a(z);
            F = P(z,z) + H(d_index,d_index);
        end
        sig=sqrt(diag(F));
        if any(diag(F)<kalman_tol) || rcond(F./(sig*sig')) < kalman_tol
            if ~all(abs(F(:))<kalman_tol)
                %use diffuse filter
                return
            else
                %pathological case, discard draw
                return
            end
        else
            F_singular = 0;
            log_dF      = log(det(F./(sig*sig')))+2*sum(log(sig));
            iF      = inv(F./(sig*sig'))./(sig*sig');
            lik(s) = log_dF + transpose(v)*iF*v + length(d_index)*log(2*pi);
            if Zflag
                K = P*z'*iF;
                P = T*(P-K*z*P)*transpose(T)+QQ;
            else
                K = P(:,z)*iF;
                P = T*(P-K*P(z,:))*transpose(T)+QQ;
            end
            a = T*(a+K*v);
            if t>=no_more_missing_observations
                notsteady = max(abs(K(:)-oldK))>riccati_tol;
                oldK = K(:);
            end
        end
    end
    t = t+1;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Divide by two.
lik(1:s) = .5*lik(1:s);

% Call steady state Kalman filter if needed.
if t<=last
    [tmp, lik(s+1:end)] = kalman_filter_ss(Y,t,last,a,T,K,iF,log_dF,Z,pp,Zflag);
end

% Compute minus the log-likelihood.
if presample>=diffuse_periods
    LIK = sum(lik(1+presample-diffuse_periods:end));
else
    LIK = sum(lik);
end