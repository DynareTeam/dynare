function [dLIK,dlik,a,Pstar] = kalman_filter_d(Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, riccati_tol, presample, T, R, Q, H, Z, mm, pp, rr)
% Computes the diffuse likelihood of a state space model.
%
% INPUTS
%    Y           [double]      pp*smpl matrix of (detrended) data, where pp is the number of observed variables.
%    start       [integer]     scalar, first observation.
%    last        [integer]     scalar, last observation.
%    a           [double]      mm*1 vector, levels of the state variables.
%    Pinf        [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    Pstar       [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    kalman_tol  [double]      scalar, tolerance parameter (rcond) of F_star.
%    diffuse_kalman_tol [double]      scalar, tolerance parameter (rcond) of Pinf to signify end of diffuse filtering and Finf.
%    riccati_tol [double]      scalar, tolerance parameter (riccati iteration);
%                              not used in this filter as usually diffuse phase will be left before convergence of filter to steady state.
%    presample   [integer]     scalar, presampling if strictly positive.
%    T           [double]      mm*mm matrix, transition matrix in  the state equations.
%    R           [double]      mm*rr matrix relating the structural innovations to the state vector.
%    Q           [double]      rr*rr covariance matrix of the structural innovations.
%    H           [double]      pp*pp covariance matrix of the measurement errors (if H is equal to zero (scalar) there is no measurement error).
%    Z           [double]      pp*mm matrix, selection matrix or pp linear independent combinations of the state vector.
%    mm          [integer]     scalar, number of state variables.
%    pp          [integer]     scalar, number of observed variables.
%    rr          [integer]     scalar, number of structural innovations.
%
% OUTPUTS
%    LIK         [double]      scalar, minus loglikelihood
%    lik         [double]      smpl*1 vector, log density of each vector of observations.
%    a           [double]      mm*1 vector, current estimate of the state vector.
%    Pstar       [double]      mm*mm matrix, covariance matrix of the state vector.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003), in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98.
%   and
%   Durbin/Koopman (2012): "Time Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 5 and 7.2

% Copyright (C) 2004-2017 Dynare Team
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
s    = 0;

while rank(Z*Pinf*Z',diffuse_kalman_tol) && (t<=last)
    s = t-start+1;
    v = Y(:,t)-Z*a;                                                     %get prediction error v^(0) in (5.13) DK (2012)
    Finf  = Z*Pinf*Z';                                                  % (5.7) in DK (2012)
                                                                        %do case distinction based on whether F_{\infty,t} has full rank or 0 rank
    if rcond(Finf) < diffuse_kalman_tol                                 %F_{\infty,t} = 0
        if ~all(abs(Finf(:)) < diffuse_kalman_tol)                      %rank-deficient but not rank 0
                                                                        % The univariate diffuse kalman filter should be used instead.
            return
        else                                                            %rank of F_{\infty,t} is 0
            Fstar  = Z*Pstar*Z' + H;                                    % (5.7) in DK (2012)
            if rcond(Fstar) < kalman_tol                                %F_{*} is singular
                if ~all(abs(Fstar(:))<kalman_tol)
                    % The univariate diffuse kalman filter should be used.
                    return
                else                                                    %rank 0
                                                                        %pathological case, discard draw
                    return
                end
            else
                iFstar = inv(Fstar);
                dFstar = det(Fstar);
                Kstar  = Pstar*Z'*iFstar;                               %(5.15) of DK (2012) with Kstar=T^{-1}*K^(0)
                dlik(s)= log(dFstar) + v'*iFstar*v;                     %set w_t to bottom case in bottom equation page 172, DK (2012)
                Pinf   = T*Pinf*transpose(T);                           % (5.16) DK (2012)
                Pstar  = T*(Pstar-Pstar*Z'*Kstar')*T'+QQ;               % (5.17) DK (2012)
                a      = T*(a+Kstar*v);                                 % (5.13) DK (2012)
            end
        end
    else                                                                %F_{\infty,t} positive definite
                                                                        %To compare to DK (2012), this block makes use of the following transformation
                                                                        %Kstar=T^{-1}*K^{(1)}=M_{*}*F^{(1)}+M_{\infty}*F^{(2)}
                                                                        %     =P_{*}*Z'*F^{(1)}+P_{\infty}*Z'*((-1)*(-F_{\infty}^{-1})*F_{*}*(F_{\infty}^{-1}))
                                                                        %     =[P_{*}*Z'-Kinf*F_{*})]*F^{(1)}
                                                                        %Make use of L^{0}'=(T-K^{(0)}*Z)'=(T-T*M_{\infty}*F^{(1)}*Z)'
                                                                        %                  =(T-T*P_{\infty*Z'*F^{(1)}*Z)'=(T-T*Kinf*Z)'
                                                                        %                  = (T*(I-*Kinf*Z))'=(I-Z'*Kinf')*T'
                                                                        %P_{*}=T*P_{\infty}*L^{(1)}+T*P_{*}*L^{(0)}+RQR
                                                                        %     =T*[(P_{\infty}*(-K^{(1)*Z}))+P_{*}*(I-Z'*Kinf')*T'+RQR]
        dlik(s)= log(det(Finf));                                        %set w_t to top case in bottom equation page 172, DK (2012)
        iFinf  = inv(Finf);
        Kinf   = Pinf*Z'*iFinf;                                         %define Kinf=T^{-1}*K_0 with M_{\infty}=Pinf*Z'
        Fstar  = Z*Pstar*Z' + H;                                        %(5.7) DK(2012)
        Kstar  = (Pstar*Z'-Kinf*Fstar)*iFinf;                           %(5.12) DK(2012); note that there is a typo in DK (2003) with "+ Kinf" instead of "- Kinf", but it is correct in their appendix
        Pstar  = T*(Pstar-Pstar*Z'*Kinf'-Pinf*Z'*Kstar')*T'+QQ;         %(5.14) DK(2012)
        Pinf   = T*(Pinf-Pinf*Z'*Kinf')*T';                             %(5.14) DK(2012)
        a      = T*(a+Kinf*v);                                          %(5.13) DK(2012)
    end
    t = t+1;
end

if t>last
    warning(['kalman_filter_d: There isn''t enough information to estimate the initial conditions of the nonstationary variables. The diffuse Kalman filter never left the diffuse stage.']);
    dLIK = NaN;
    return
end

dlik = dlik(1:s);
dlik = .5*(dlik + pp*log(2*pi));

dLIK = sum(dlik(1+presample:end));
