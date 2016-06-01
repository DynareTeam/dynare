function [fhat,xhat,fcount,retcode] = csminit1(fcn,x0,penalty,f0,g0,badg,H0,Verbose,varargin)
% [fhat,xhat,fcount,retcode] = csminit1(fcn,x0,penalty,f0,g0,badg,H0,Verbose,varargin)
% 
% Inputs: 
%   fcn:        [string]        string naming the objective function to be minimized
%   x0:         [npar by 1]     initial value of the parameter vector
%   penalty:    [scalar]        variable penalty in case of failure of objective function
%   f0:         [scalar]        initial value of the fucntion
%   g0:         [npar by 1]     initial value of the gradient vector
%   badg        [scalar]        indicator for problem in gradient computation
%   H0:         [npar by npar]  initial value for the inverse Hessian.  Must be positive definite.
%   Verbose:    [scalar]        indicator for verbose computations
%   varargin:                   Optional additional inputs that get handed off to fcn each
%                               time it is called.

% Outputs:
%   fhat:       [scalar]        function value at minimum
%   xhat:       [npar by 1]     parameter vector at minimum
%   fcount      [scalar]        function iteration count upon termination
%   retcode     [scalar]        0: normal step
%                               1: zero gradient.
%                               5: largest step still improves too fast.
%                               2,4: back and forth adjustment of stepsize didn't finish.  
%                               3: smallest stepsize still improves too slow
%                               6: no improvement found
%---------------------
% Modified 7/22/96 to omit variable-length P list, for efficiency and compilation.
% Places where the number of P's need to be altered or the code could be returned to
% its old form are marked with ARGLIST comments.
%
% Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
% update.
%
% Fixed 7/19/93 to flip eigenvalues of H to get better performance when
% it's not psd.
% 
% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/csminit.m
% 
% Copyright (C) 1993-2007 Christopher Sims
% Copyright (C) 2008-2016 Dynare Team
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

%tailstr = ')';
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%end
%ANGLE = .03;
ANGLE = .005;
%THETA = .03;
THETA = .3; %(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
FCHANGE = 1000;
MINLAMB = 1e-9;
% fixed 7/15/94
% MINDX = .0001;
% MINDX = 1e-6;
MINDFAC = .01;
fcount=0;
lambda=1;
xhat=x0;
f=f0;
fhat=f0;
g = g0;
gnorm = norm(g);
%
if (gnorm < 1.e-12) && ~badg % put ~badg 8/4/94
    retcode =1;
    dxnorm=0;
    % gradient convergence
else
    % with badg true, we don't try to match rate of improvement to directional
    % derivative.  We're satisfied just to get some improvement in f.
    %
    %if(badg)
    %   dx = -g*FCHANGE/(gnorm*gnorm);
    %  dxnorm = norm(dx);
    %  if dxnorm > 1e12
    %     disp('Bad, small gradient problem.')
    %     dx = dx*FCHANGE/dxnorm;
    %   end
    %else
    % Gauss-Newton step;
    %---------- Start of 7/19/93 mod ---------------
    %[v d] = eig(H0);
    %toc
    %d=max(1e-10,abs(diag(d)));
    %d=abs(diag(d));
    %dx = -(v.*(ones(size(v,1),1)*d'))*(v'*g);
    %      toc
    dx = -H0*g;
    %      toc
    dxnorm = norm(dx);
    if dxnorm > 1e12
        disp_verbose('Near-singular H problem.',Verbose)
        dx = dx*FCHANGE/dxnorm;
    end
    dfhat = dx'*g0;
    %end
    %
    %
    if ~badg
        % test for alignment of dx with gradient and fix if necessary
        a = -dfhat/(gnorm*dxnorm);
        if a<ANGLE
            dx = dx - (ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g;
            dfhat = dx'*g;
            dxnorm = norm(dx);
            disp_verbose(sprintf('Correct for low angle: %g',a),Verbose)
        end
    end
        disp_verbose(sprintf('Predicted improvement: %18.9f',-dfhat/2),Verbose)
    %
    % Have OK dx, now adjust length of step (lambda) until min and
    % max improvement rate criteria are met.
    done=0;
    factor=3;
    shrink=1;
    lambdaMin=0;
    lambdaMax=inf;
    lambdaPeak=0;
    fPeak=f0;
    lambdahat=0;
    while ~done
        if size(x0,2)>1
            dxtest=x0+dx'*lambda;
        else
            dxtest=x0+dx*lambda;
        end
        % home
        f = penalty_objective_function(dxtest,fcn,penalty,varargin{:});
        %ARGLIST
        %f = feval(fcn,dxtest,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
        % f = feval(fcn,x0+dx*lambda,P1,P2,P3,P4,P5,P6,P7,P8);
            disp_verbose(sprintf('lambda = %10.5g; f = %20.7f',lambda,f ),Verbose)
        %debug
        %disp(sprintf('Improvement too great? f0-f: %g, criterion: %g',f0-f,-(1-THETA)*dfhat*lambda))
        if f<fhat
            fhat=f;
            xhat=dxtest;
            lambdahat = lambda;
        end
        fcount=fcount+1;
        shrinkSignal = (~badg & (f0-f < max([-THETA*dfhat*lambda 0]))) | (badg & (f0-f) < 0) ;
        growSignal = ~badg & ( (lambda > 0)  &  (f0-f > -(1-THETA)*dfhat*lambda) );
        if  shrinkSignal  &&   ( (lambda>lambdaPeak) || (lambda<0) )
            if (lambda>0) && ((~shrink) || (lambda/factor <= lambdaPeak))
                shrink=1;
                factor=factor^.6;
                while lambda/factor <= lambdaPeak
                    factor=factor^.6;
                end
                %if (abs(lambda)*(factor-1)*dxnorm < MINDX) || (abs(lambda)*(factor-1) < MINLAMB)
                if abs(factor-1)<MINDFAC
                    if abs(lambda)<4
                        retcode=2;
                    else
                        retcode=7;
                    end
                    done=1;
                end
            end
            if (lambda<lambdaMax) && (lambda>lambdaPeak)
                lambdaMax=lambda;
            end
            lambda=lambda/factor;
            if abs(lambda) < MINLAMB
                if (lambda > 0) && (f0 <= fhat)
                    % try going against gradient, which may be inaccurate
                    if Verbose
                        lambda = -lambda*factor^6
                    else
                        lambda = -lambda*factor^6;                        
                    end
                else
                    if lambda < 0
                        retcode = 6;
                    else
                        retcode = 3;
                    end
                    done = 1;
                end
            end
        elseif  (growSignal && lambda>0) ||  (shrinkSignal && ((lambda <= lambdaPeak) && (lambda>0)))
            if shrink
                shrink=0;
                factor = factor^.6;
                %if ( abs(lambda)*(factor-1)*dxnorm< MINDX ) || ( abs(lambda)*(factor-1)< MINLAMB)
                if abs(factor-1)<MINDFAC
                    if abs(lambda)<4
                        retcode=4;
                    else
                        retcode=7;
                    end
                    done=1;
                end
            end
            if ( f<fPeak ) && (lambda>0)
                fPeak=f;
                lambdaPeak=lambda;
                if lambdaMax<=lambdaPeak
                    lambdaMax=lambdaPeak*factor*factor;
                end
            end
            lambda=lambda*factor;
            if abs(lambda) > 1e20;
                retcode = 5;
                done =1;
            end
        else
            done=1;
            if factor < 1.2
                retcode=7;
            else
                retcode=0;
            end
        end
    end
end

disp_verbose(sprintf('Norm of dx %10.5g', dxnorm),Verbose)
