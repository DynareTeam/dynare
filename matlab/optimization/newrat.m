function [xparam1, hh, gg, fval, igg, hess_info] = newrat(func0, x, bounds, analytic_derivation, ftol0, nit, flagg, Verbose, Save_files, hess_info, varargin)
%  [xparam1, hh, gg, fval, igg, hess_info] = newrat(func0, x, bounds, analytic_derivation, ftol0, nit, flagg, Verbose, Save_files, hess_info, varargin)
%
%  Optimiser with outer product gradient and with sequences of univariate steps
%  uses Chris Sims subroutine for line search
%
%  Inputs:
%  - func0                  name of the function that also outputs the single contributions at times t=1,...,T
%                           of the log-likelihood to compute outer product gradient
%  - x                      starting guess
%  - analytic_derivation    1 if analytic derivatives, 0 otherwise
%  - ftol0                  termination criterion for function change
%  - nit                    maximum number of iterations
%  - flagg                  Indicator how to compute final Hessian (In each iteration, Hessian is computed with outer product gradient)  
%                           0: final Hessian computed with outer product gradient
%                           1: final 'mixed' Hessian: diagonal elements computed with 
%                               numerical second order derivatives with correlation structure 
%                               as from outer product gradient
%                           2: full numerical Hessian
%  - Verbose                1 if explicit output is requested
%  - Save_files             1 if intermediate output is to be saved 
%  - hess_info              structure storing the step sizes for
%                           computation of Hessian
%  - varargin               other inputs:
%                           varargin{1} --> DynareDataset
%                           varargin{2} --> DatasetInfo
%                           varargin{3} --> DynareOptions
%                           varargin{4} --> Model
%                           varargin{5} --> EstimatedParameters
%                           varargin{6} --> BayesInfo
%                           varargin{7} --> Bounds
%                           varargin{8} --> DynareResults
% 
% Outputs
% - xparam1                 parameter vector at optimum
% - hh                      hessian
% - gg                      gradient
% - fval                    function value
% - igg                     inverted outer product hessian
% - hess_info               structure with updated step length

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

% initialize variable penalty    
penalty = 1e8;

icount=0;
nx=length(x);
xparam1=x;
%ftol0=1.e-6;
htol_base = max(1.e-7, ftol0);
flagit=0;  % mode of computation of hessian in each iteration
ftol=ftol0;
gtol=1.e-3;
htol=htol_base;
htol0=htol_base;
gibbstol=length(varargin{6}.pshape)/50; %25;

% force fcn, grad to function handle
if ischar(func0)
    func0 = str2func(func0);
end

% func0 = str2func([func2str(func0),'_hh']);
% func0 = func0;
[fval0,exit_flag,gg,hh]=penalty_objective_function(x,func0,penalty,varargin{:});
fval=fval0;

% initialize mr_gstep and mr_hessian

outer_product_gradient=1;
if isempty(hh)
    [dum, gg, htol0, igg, hhg, h1, hess_info]=mr_hessian(x,func0,penalty,flagit,htol,hess_info,varargin{:});
    if isempty(dum),
        outer_product_gradient=0;
        igg = 1e-4*eye(nx);
    else
        hh0 = reshape(dum,nx,nx);
        hh=hhg;
        if min(eig(hh0))<0
            hh0=hhg; %generalized_cholesky(hh0);
        elseif flagit==2
            hh=hh0;
            igg=inv(hh);
        end
    end
    if max(htol0)>htol
        skipline()
        disp_verbose('Numerical noise in the likelihood',Verbose)
        disp_verbose('Tolerance has to be relaxed',Verbose)
        skipline()
    end
else
    hh0=hh;
    hhg=hh;
    igg=inv(hh);
    h1=[];
end
H = igg;
disp_verbose(['Gradient norm ',num2str(norm(gg))],Verbose)
ee=eig(hh);
disp_verbose(['Minimum Hessian eigenvalue ',num2str(min(ee))],Verbose)
disp_verbose(['Maximum Hessian eigenvalue ',num2str(max(ee))],Verbose)
g=gg;
check=0;
if max(eig(hh))<0
    disp_verbose('Negative definite Hessian! Local maximum!',Verbose)
    pause
end
if Save_files
    save('m1.mat','x','hh','g','hhg','igg','fval0')
end

igrad=1;
igibbs=1;
inx=eye(nx);
jit=0;
nig=[];
ig=ones(nx,1);
ggx=zeros(nx,1);
while norm(gg)>gtol && check==0 && jit<nit
    jit=jit+1;
    tic1 = tic;
    icount=icount+1;
    penalty = fval0(icount);
    disp_verbose([' '],Verbose)
    disp_verbose(['Iteration ',num2str(icount)],Verbose)
    [fval,x0,fc,retcode] = csminit1(func0,xparam1,penalty,fval0(icount),gg,0,H,Verbose,varargin{:});
    if igrad
        [fval1,x01,fc,retcode1] = csminit1(func0,x0,penalty,fval,gg,0,inx,Verbose,varargin{:});
        if (fval-fval1)>1
            disp_verbose('Gradient step!!',Verbose)
        else
            igrad=0;
        end
        fval=fval1;
        x0=x01;
    end
    if length(find(ig))<nx
        ggx=ggx*0;
        ggx(find(ig))=gg(find(ig));
        if analytic_derivation,
            hhx=hh;
        else
            hhx = reshape(dum,nx,nx);
        end
        iggx=eye(length(gg));
        iggx(find(ig),find(ig)) = inv( hhx(find(ig),find(ig)) );
        [fvala,x0,fc,retcode] = csminit1(func0,x0,penalty,fval,ggx,0,iggx,Verbose,varargin{:});
    end
    x0 = check_bounds(x0,bounds);
    [fvala, x0, ig] = mr_gstep(h1,x0,bounds,func0,penalty,htol0,Verbose,Save_files,varargin{:});
    x0 = check_bounds(x0,bounds);
    nig=[nig ig];
    disp_verbose('Sequence of univariate steps!!',Verbose)
    fval=fvala;
    if (fval0(icount)-fval)<ftol && flagit==0
        disp_verbose('Try diagonal Hessian',Verbose)
        ihh=diag(1./(diag(hhg)));
        [fval2,x0,fc,retcode2] = csminit1(func0,x0,penalty,fval,gg,0,ihh,Verbose,varargin{:});
        x0 = check_bounds(x0,bounds);
        if (fval-fval2)>=ftol
            disp_verbose('Diagonal Hessian successful',Verbose)
        end
        fval=fval2;
    end
    if (fval0(icount)-fval)<ftol && flagit==0
        disp_verbose('Try gradient direction',Verbose)
        ihh0=inx.*1.e-4;
        [fval3,x0,fc,retcode3] = csminit1(func0,x0,penalty,fval,gg,0,ihh0,Verbose,varargin{:});
        x0 = check_bounds(x0,bounds);
        if (fval-fval3)>=ftol
            disp_verbose('Gradient direction successful',Verbose)
        end
        fval=fval3;
    end
    xparam1=x0;
    x(:,icount+1)=xparam1;
    fval0(icount+1)=fval;
    if (fval0(icount)-fval)<ftol
        disp_verbose('No further improvement is possible!',Verbose)
        check=1;
        if analytic_derivation,
            [fvalx,exit_flag,gg,hh]=penalty_objective_function(xparam1,func0,penalty,varargin{:});
            hhg=hh;
            H = inv(hh);
        else
        if flagit==2
            hh=hh0;
        elseif flagg>0
            [dum, gg, htol0, igg, hhg, h1, hess_info]=mr_hessian(xparam1,func0,penalty,flagg,ftol0,hess_info,varargin{:});
            if flagg==2
                hh = reshape(dum,nx,nx);
                ee=eig(hh);
                if min(ee)<0
                    hh=hhg;
                end
            else
                hh=hhg;
            end
        end
        end
        disp_verbose(['Actual dxnorm ',num2str(norm(x(:,end)-x(:,end-1)))],Verbose)
        disp_verbose(['FVAL          ',num2str(fval)],Verbose)
        disp_verbose(['Improvement   ',num2str(fval0(icount)-fval)],Verbose)
        disp_verbose(['Ftol          ',num2str(ftol)],Verbose)
        disp_verbose(['Htol          ',num2str(max(htol0))],Verbose)
        disp_verbose(['Gradient norm  ',num2str(norm(gg))],Verbose)
        ee=eig(hh);
        disp_verbose(['Minimum Hessian eigenvalue ',num2str(min(ee))],Verbose)
        disp_verbose(['Maximum Hessian eigenvalue ',num2str(max(ee))],Verbose)
        g(:,icount+1)=gg;
    else
        df = fval0(icount)-fval;
        disp_verbose(['Actual dxnorm ',num2str(norm(x(:,end)-x(:,end-1)))],Verbose)
        disp_verbose(['FVAL          ',num2str(fval)],Verbose)
        disp_verbose(['Improvement   ',num2str(df)],Verbose)
        disp_verbose(['Ftol          ',num2str(ftol)],Verbose)
        disp_verbose(['Htol          ',num2str(max(htol0))],Verbose)
        htol=htol_base;
        if norm(x(:,icount)-xparam1)>1.e-12 && analytic_derivation==0,
            try
                if Save_files
                    save('m1.mat','x','fval0','nig','-append')
                end
            catch
                if Save_files
                    save('m1.mat','x','fval0','nig')
                end
            end
            [dum, gg, htol0, igg, hhg, h1, hess_info]=mr_hessian(xparam1,func0,penalty,flagit,htol,hess_info,varargin{:});
            if isempty(dum),
                outer_product_gradient=0;
            end
            if max(htol0)>htol
                skipline()
                disp_verbose('Numerical noise in the likelihood',Verbose)
                disp_verbose('Tolerance has to be relaxed',Verbose)
                skipline()
            end
            if ~outer_product_gradient,
                H = bfgsi1(H,gg-g(:,icount),xparam1-x(:,icount),Verbose,Save_files);
                hh=inv(H);
                hhg=hh;
            else
                hh0 = reshape(dum,nx,nx);
                hh=hhg;
                if flagit==2
                    if min(eig(hh0))<=0
                        hh0=hhg; %generalized_cholesky(hh0);
                    else
                        hh=hh0;
                        igg=inv(hh);
                    end
                end
                H = igg;
            end
        elseif analytic_derivation,
            [fvalx,exit_flag,gg,hh]=penalty_objective_function(xparam1,func0,penalty,varargin{:});
            hhg=hh;
            H = inv(hh);
        end
        disp_verbose(['Gradient norm  ',num2str(norm(gg))],Verbose)
        ee=eig(hh);
        disp_verbose(['Minimum Hessian eigenvalue ',num2str(min(ee))],Verbose)
        disp_verbose(['Maximum Hessian eigenvalue ',num2str(max(ee))],Verbose)
        if max(eig(hh))<0
            disp_verbose('Negative definite Hessian! Local maximum!',Verbose)
            pause(1)
        end
        t=toc(tic1);
        disp_verbose(['Elapsed time for iteration ',num2str(t),' s.'],Verbose)
        g(:,icount+1)=gg;
        if Save_files
            save('m1.mat','x','hh','g','hhg','igg','fval0','nig','H')
        end
    end
end
if Save_files
    save('m1.mat','x','hh','g','hhg','igg','fval0','nig')
end
if ftol>ftol0
    skipline()
    disp_verbose('Numerical noise in the likelihood',Verbose)
    disp_verbose('Tolerance had to be relaxed',Verbose)
    skipline()
end

if jit==nit
    skipline()
    disp_verbose('Maximum number of iterations reached',Verbose)
    skipline()
end

if norm(gg)<=gtol
    disp_verbose(['Estimation ended:'],Verbose)
    disp_verbose(['Gradient norm < ', num2str(gtol)],Verbose)
end
if check==1,
    disp_verbose(['Estimation successful.'],Verbose)
end

return


function x = check_bounds(x,bounds)

inx = find(x>=bounds(:,2));
if ~isempty(inx),
    x(inx) = bounds(inx,2)-eps;
end

inx = find(x<=bounds(:,1));
if ~isempty(inx),
    x(inx) = bounds(inx,1)+eps;
end
