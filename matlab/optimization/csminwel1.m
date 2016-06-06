function [fh,xh,gh,H,itct,fcount,retcodeh] = csminwel1(fcn,x0,H0,grad,crit,nit,method,epsilon,Verbose,Save_files,varargin)
%[fhat,xhat,ghat,Hhat,itct,fcount,retcodeh] = csminwel1(fcn,x0,H0,grad,crit,nit,method,epsilon,varargin)
% Inputs:
%   fcn:    [string]        string naming the objective function to be minimized
%   x0:     [npar by 1]     initial value of the parameter vector
%   H0:     [npar by npar]  initial value for the inverse Hessian.  Must be positive definite.
%   grad:   [string or empty matrix] Either a string naming a function that calculates the gradient, or the null matrix.
%                           If it's null, the program calculates a numerical gradient.  In this case fcn must
%                           be written so that it can take a matrix argument and produce a row vector of values.
%   crit:   [scalar]        Convergence criterion.  Iteration will cease when it proves impossible to improve the
%                           function value by more than crit.
%   nit:    [scalar]        Maximum number of iterations.
%   method: [scalar]        integer scalar for selecting gradient method: 2, 3 or 5 points formula.
%   epsilon: [scalar]       scalar double, numerical differentiation increment
%   varargin:               Optional additional inputs that get handed off to fcn each
%                           time it is called.
% 
%        Note that if the program ends abnormally, it is possible to retrieve the current x,
%        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
%        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
%        writes g2.mat and g3.mat as well. If all were written at about the same time, any of them
%        may be a decent starting point. One can also start from the one with best function value.)
% 
% Outputs:
%   fh:     [scalar]        function value at minimum
%   xh:     [npar by 1]     parameter vector at minimum
%   gh      [npar by 1]     gradient vector
%   H       [npar by npar]  inverse of the Hessian matrix
%   itct    [scalar]        iteration count upon termination
%   fcount  [scalar]        function iteration count upon termination
%   retcodeh [scalar]       return code:
%                               0: normal step
%                               1: zero gradient 
%                               2: back and forth on step length never finished
%                               3: smallest step still improving too slow
%                               4: back and forth on step length never finished
%                               5: largest step still improving too fast
%                               6: smallest step still improving too slow, reversed gradient
%                               7: warning: possible inaccuracy in H matrix
% 
% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/csminwel.m
% 
% Copyright (C) 1993-2007 Christopher Sims
% Copyright (C) 2006-2015 Dynare Team
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

% initialize variable penalty    
penalty = 1e8;
fh = [];
xh = [];
[nx,no]=size(x0);
nx=max(nx,no);
NumGrad= isempty(grad);
done=0;
itct=0;
fcount=0;
gh = [];
H = [];
retcodeh = [];

% force fcn, grad to function handle
if ischar(fcn)
    fcn = str2func(fcn);
end
if ischar(grad)
    grad = str2func(grad);
end
%tailstr = ')';
%stailstr = [];
% Lines below make the number of Pi's optional.  This is inefficient, though, and precludes
% use of the matlab compiler.  Without them, we use feval and the number of Pi's must be
% changed with the editor for each application.  Places where this is required are marked
% with ARGLIST comments
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%   stailstr=[' P' num2str(i) stailstr];
%end

[f0,cost_flag,arg1] = penalty_objective_function(x0,fcn,penalty,varargin{:});

if ~cost_flag
    disp_verbose('Bad initial parameter.',Verbose)
    return
end

if NumGrad
    [g, badg]=get_num_grad(method,fcn,penalty,f0,x0,epsilon,varargin{:});
elseif ischar(grad)
    [g,badg] = grad(x0,varargin{:});
else
    g=arg1;
    badg=0;
end
retcode3=101;
x=x0;
f=f0;
H=H0;
cliff=0;
while ~done
    % penalty for dsge_likelihood and dsge_var_likelihood
    penalty = f;

    g1=[]; g2=[]; g3=[];
    %addition fj. 7/6/94 for control
    disp_verbose('-----------------',Verbose)
    disp_verbose(sprintf('f at the beginning of new iteration, %20.10f',f),Verbose)
    %-----------Comment out this line if the x vector is long----------------
    %   disp_verbose([sprintf('x = ') sprintf('%15.8g %15.8g %15.8g %15.8g\n',x)]);
    %-------------------------
    itct=itct+1;
    [f1, x1, fc, retcode1] = csminit1(fcn,x,penalty,f,g,badg,H,Verbose,varargin{:});
    fcount = fcount+fc;
    % erased on 8/4/94
    % if (retcode == 1) || (abs(f1-f) < crit)
    %    done=1;
    % end
    % if itct > nit
    %    done = 1;
    %    retcode = -retcode;
    % end
    if retcode1 ~= 1
        if retcode1==2 || retcode1==4
            wall1=1; badg1=1;
        else
            if NumGrad
                [g1, badg1]=get_num_grad(method,fcn,penalty,f1,x1,epsilon,varargin{:});
            elseif ischar(grad),
                [g1, badg1] = grad(x1,varargin{:});
            else
                [junk1,cost_flag,g1] = penalty_objective_function(x1,fcn,penalty,varargin{:});
                badg1 = ~cost_flag;
            end
            wall1=badg1;
            % g1
            if Save_files
                save('g1.mat','g1','x1','f1','varargin');
            end
        end
        if wall1 % && (~done) by Jinill
                 % Bad gradient or back and forth on step length.  Possibly at
                 % cliff edge.  Try perturbing search direction.
                 %
                 %fcliff=fh;xcliff=xh;
            Hcliff=H+diag(diag(H).*rand(nx,1));
            disp_verbose('Cliff.  Perturbing search direction.',Verbose)
            [f2, x2, fc, retcode2] = csminit1(fcn,x,penalty,f,g,badg,Hcliff,Verbose,varargin{:});
            fcount = fcount+fc; % put by Jinill
            if  f2 < f
                if retcode2==2 || retcode2==4
                    wall2=1; badg2=1;
                else
                    if NumGrad
                        [g2, badg2]=get_num_grad(method,fcn,penalty,f2,x2,epsilon,varargin{:});
                    elseif ischar(grad),
                        [g2, badg2] = grad(x2,varargin{:});
                    else
                        [junk2,cost_flag,g2] = penalty_objective_function(x1,fcn,penalty,varargin{:});
                        badg2 = ~cost_flag;
                    end
                    wall2=badg2;
                    % g2
                    if Verbose
                        badg2
                    end
                    if Save_files
                        save('g2.mat','g2','x2','f2','varargin');
                    end
                end
                if wall2
                    disp_verbose('Cliff again.  Try traversing',Verbose)
                    if norm(x2-x1) < 1e-13
                        f3=f; x3=x; badg3=1;retcode3=101;
                    else
                        gcliff=((f2-f1)/((norm(x2-x1))^2))*(x2-x1);
                        if(size(x0,2)>1)
                            gcliff=gcliff';
                        end
                        [f3, x3, fc, retcode3] = csminit1(fcn,x,penalty,f,gcliff,0,eye(nx),Verbose,varargin{:});
                        fcount = fcount+fc; % put by Jinill
                        if retcode3==2 || retcode3==4
                            wall3=1; 
                            badg3=1;
                        else
                            if NumGrad
                                [g3, badg3]=get_num_grad(method,fcn,penalty,f3,x3,epsilon,varargin{:});
                            elseif ischar(grad),
                                [g3, badg3] = grad(x3,varargin{:});
                            else
                                [junk3,cost_flag,g3] = penalty_objective_function(x1,fcn,penalty,varargin{:});
                                badg3 = ~cost_flag;
                            end
                            wall3=badg3;
                            % g3
                            if Save_files
                                save('g3.mat','g3','x3','f3','varargin');
                            end
                        end
                    end
                else
                    f3=f; x3=x; badg3=1; retcode3=101;
                end
            else
                f3=f; x3=x; badg3=1;retcode3=101;
            end
        else
            % normal iteration, no walls, or else we're finished here.
            f2=f; f3=f; badg2=1; badg3=1; retcode2=101; retcode3=101;
        end
    else
        f2=f;f3=f;f1=f;retcode2=retcode1;retcode3=retcode1;
    end
    %how to pick gh and xh
    if f3 < f - crit && badg3==0 && f3 < f2 && f3 < f1
        ih=3;
        fh=f3;xh=x3;gh=g3;badgh=badg3;retcodeh=retcode3;
    elseif f2 < f - crit && badg2==0 && f2 < f1
        ih=2;
        fh=f2;xh=x2;gh=g2;badgh=badg2;retcodeh=retcode2;
    elseif f1 < f - crit && badg1==0
        ih=1;
        fh=f1;xh=x1;gh=g1;badgh=badg1;retcodeh=retcode1;
    else
        [fh,ih] = min([f1,f2,f3]);
        %disp_verbose(sprintf('ih = %d',ih))
        %eval(['xh=x' num2str(ih) ';'])
        switch ih
          case 1
            xh=x1;
          case 2
            xh=x2;
          case 3
            xh=x3;
        end %case
            %eval(['gh=g' num2str(ih) ';'])
            %eval(['retcodeh=retcode' num2str(ih) ';'])
        retcodei=[retcode1,retcode2,retcode3];
        retcodeh=retcodei(ih);
        if exist('gh')
            nogh=isempty(gh);
        else
            nogh=1;
        end
        if nogh
            if NumGrad
                [gh, badgh]=get_num_grad(method,fcn,penalty,fh,xh,epsilon,varargin{:});
            elseif ischar(grad),
                [gh, badgh] = grad(xh,varargin{:});
            else
                [junkh,cost_flag,gh] = penalty_objective_function(x1,fcn,penalty,varargin{:});
                badgh = ~cost_flag;
            end
        end
        badgh=1;
    end
    %end of picking
    stuck = (abs(fh-f) < crit);
    if (~badg) && (~badgh) && (~stuck)
        H = bfgsi1(H,gh-g,xh-x,Verbose,Save_files);
    end
    disp_verbose('----',Verbose)
    disp_verbose(sprintf('Improvement on iteration %d = %18.9f',itct,f-fh),Verbose)
    % if Verbose
    if itct > nit
        disp_verbose('iteration count termination',Verbose)
        done = 1;
    elseif stuck
        disp_verbose('improvement < crit termination',Verbose)
        done = 1;
    end
    rc=retcodeh;
    if Verbose || done
        if rc ==0
            %do nothing, just a normal step
        elseif rc == 1
            disp_verbose('zero gradient',Verbose)
        elseif rc == 6
            disp_verbose('smallest step still improving too slow, reversed gradient',Verbose)
        elseif rc == 5
            disp_verbose('largest step still improving too fast',Verbose)
        elseif (rc == 4) || (rc==2)
            disp_verbose('back and forth on step length never finished',Verbose)
        elseif rc == 3
            disp_verbose('smallest step still improving too slow',Verbose)
        elseif rc == 7
            disp_verbose('warning: possible inaccuracy in H matrix',Verbose)
        else
            error('Unaccounted Case, please contact the developers',Verbose)
        end
    end
     
    f=fh;
    x=xh;
    g=gh;
    badg=badgh;
end

end

function [g, badg]=get_num_grad(method,fcn,penalty,f0,x0,epsilon,varargin)
    switch method
      case 2
        [g,badg] = numgrad2(fcn, f0, x0, penalty, epsilon, varargin{:});
      case 3
        [g,badg] = numgrad3(fcn, f0, x0, penalty, epsilon, varargin{:});
      case 5
        [g,badg] = numgrad5(fcn, f0, x0, penalty, epsilon, varargin{:});
      case 13
        [g,badg] = numgrad3_(fcn, f0, x0, penalty, epsilon, varargin{:});
      case 15
        [g,badg] = numgrad5_(fcn, f0, x0, penalty, epsilon, varargin{:});
      otherwise
        error('csminwel1: Unknown method for gradient evaluation!')
    end
end