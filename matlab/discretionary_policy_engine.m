function [H,G,retcode]=discretionary_policy_engine(AAlag,AA0,AAlead,BB,bigw,instr_id,beta,solve_maxit,discretion_tol,qz_criterium,H00,verbose)

% Solves the discretionary problem for a model of the form:
%
%   Loss=E_0 sum_{t=0}^{\infty} beta^t [y_t'*W*y+x_t'*Q*x_t]
%   subject to
%   AAlag*yy_{t-1}+AA0*yy_t+AAlead*yy_{t+1}+BB*e=0
%
% with W the weight on the variables in vector y_t.
%
% The solution takes the form
%   y_t=H*y_{t-1}+G*e_t
% where H=[H1;F1] and G=[H2;F2].
%
% We use the Dennis (2007, Macroeconomic Dynamics) algorithm and so we need
% to re-write the model in the form
%  A0*y_t=A1*y_{t-1}+A2*y_{t+1}+A3*x_t+A4*x_{t+1}+A5*e_t, with W the
% weight on the y_t vector and Q the weight on the x_t vector of
% instruments.
%
% Inputs:
%   AAlag               [double]    matrix of coefficients on lagged
%                                   variables
%   AA0                 [double]    matrix of coefficients on
%                                   contemporaneous variables
%   AAlead              [double]    matrix of coefficients on
%                                   leaded variables
%   BB                  [double]    matrix of coefficients on
%                                   shocks
%   bigw                [double]    matrix of coefficients on variables in
%                                   loss/objective function; stacks [W and Q]
%   instr_id            [double]    location vector of the instruments in the yy_t vector.
%   beta                [scalar]    planner discount factor
%   solve_maxit         [scalar]    maximum number of iterations
%   discretion_tol      [scalar]    convergence criterion for solution
%   qz_criterium        [scalar]    tolerance for QZ decomposition
%   H00
%   verbose             [scalar]    dummy to control verbosity
%
% Outputs:
%   H                   [double]    (endo_nbr*endo_nbr) solution matrix for endogenous
%                                   variables, stacks [H1 and H1]
%   G                   [double]    (endo_nbr*exo_nbr) solution matrix for shocks, stacks [H2 and F2]
%
%   retcode             [scalar]    return code
%
% Algorithm:
%  Dennis, Richard (2007): Optimal policy in rational expectations models: new solution algorithms,
%       Macroeconomic Dynamics, 11, 3155.

% Copyright (C) 2007-2017 Dynare Team
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

if nargin<12
    verbose=0;
    if nargin<11
        H00=[];
        if nargin<10
            qz_criterium=1.000001;
            if nargin<9
                discretion_tol=sqrt(eps);
                if nargin<8
                    solve_maxit=3000;
                    if nargin<7
                        beta=.99;
                        if nargin<6
                            error([mfilename,':: Insufficient number of input arguments'])
                        elseif nargin>12
                            error([mfilename,':: Number of input arguments cannot exceed 12'])
                        end
                    end
                end
            end
        end
    end
end

[A0,A1,A2,A3,A4,A5,W,Q,endo_nbr,exo_nbr,aux,endo_augm_id]=GetDennisMatrices(AAlag,AA0,AAlead,BB,bigw,instr_id);
% aux is a logical index of the instruments which appear with lags in the
% model. Their location in the state vector is instr_id(aux);
% endo_augm_id is index (not logical) of locations of the augmented vector
% of non-instrumental variables

AuxiliaryVariables_nbr=sum(aux);
H0=zeros(endo_nbr+AuxiliaryVariables_nbr);
if ~isempty(H00)
    H0(1:endo_nbr,1:endo_nbr)=H00;
    clear H00
end

H10=H0(endo_augm_id,endo_augm_id);
F10=H0(instr_id,endo_augm_id);

iter=0;
H1=H10;
F1=F10;
%solve equations (20) and (22) via fixed point iteration
while 1
    iter=iter+1;
    P=SylvesterDoubling(W+beta*F1'*Q*F1,beta*H1',H1,discretion_tol,solve_maxit);
    if any(any(isnan(P))) || any(any(isinf(P))) 
        P=SylvesterHessenbergSchur(W+beta*F1'*Q*F1,beta*H1',H1);
        if any(any(isnan(P))) || any(any(isinf(P)))
            retcode=2;
            return
        end
    end
    D=A0-A2*H1-A4*F1; %equation (20)
    Dinv=inv(D);
    A3DPD=A3'*Dinv'*P*Dinv;
    F1=-(Q+A3DPD*A3)\(A3DPD*A1); %component of (26)
    H1=Dinv*(A1+A3*F1); %component of (27)

    [rcode,NQ]=CheckConvergence([H1;F1]-[H10;F10],iter,solve_maxit,discretion_tol);
    if rcode
        break
    else
        if verbose
            disp(NQ)
        end
    end
    H10=H1;
    F10=F1;
end

%check if successful
retcode = 0;
switch rcode
  case 3 % nan
    retcode=63;
    retcode(2)=10000;
    if verbose
        disp([mfilename,':: NaN elements in the solution'])
    end
  case 2% maxiter
    retcode = 61;
    if verbose
        disp([mfilename,':: Maximum Number of Iterations reached'])
    end
  case 1
    BadEig=max(abs(eig(H1)))>qz_criterium;
    if BadEig
        retcode=62;
        retcode(2)=100*max(abs(eig(H1)));
        if verbose
            disp([mfilename,':: Some eigenvalues greater than qz_criterium, Model potentially unstable'])
        end
    end
end

if retcode(1)
    H=[];
    G=[];
else
    F2=-(Q+A3DPD*A3)\(A3DPD*A5); %equation (29)
    H2=Dinv*(A5+A3*F2); %equation (31)
    H=zeros(endo_nbr+AuxiliaryVariables_nbr);
    G=zeros(endo_nbr+AuxiliaryVariables_nbr,exo_nbr);
    H(endo_augm_id,endo_augm_id)=H1;
    H(instr_id,endo_augm_id)=F1;
    G(endo_augm_id,:)=H2;
    G(instr_id,:)=F2;

    % Account for auxilliary variables
    H(:,instr_id(aux))=H(:,end-(AuxiliaryVariables_nbr-1:-1:0));
    H=H(1:endo_nbr,1:endo_nbr);
    G=G(1:endo_nbr,:);
end

end


function [rcode,NQ]=CheckConvergence(Q,iter,MaxIter,crit)

NQ=max(max(abs(Q)));% norm(Q); seems too costly
if isnan(NQ)
    rcode=3;
elseif iter>MaxIter;
    rcode=2;
elseif NQ<crit
    rcode=1;
else
    rcode=0;
end

end

function [A00,A11,A22,A33,A44,A55,WW,Q,endo_nbr,exo_nbr,aux,endo_augm_id]=GetDennisMatrices(AAlag,AA0,AAlead,BB,bigw,instr_id)
%get the matrices to use the Dennis (2007) algorithm
[eq_nbr,endo_nbr]=size(AAlag);
exo_nbr=size(BB,2);
y=setdiff(1:endo_nbr,instr_id);
instr_nbr=numel(instr_id);

A0=AA0(:,y);
A1=-AAlag(:,y);
A2=-AAlead(:,y);
A3=-AA0(:,instr_id);
A4=-AAlead(:,instr_id);
A5=-BB;
W=bigw(y,y);
Q=bigw(instr_id,instr_id);
% Adjust for possible lags in instruments by creating auxiliary equations
A6=-AAlag(:,instr_id);
aux=any(A6);
AuxiliaryVariables_nbr=sum(aux);
ny=eq_nbr;
m=eq_nbr+AuxiliaryVariables_nbr;
A00=zeros(m);A00(1:ny,1:ny)=A0;A00(ny+1:end,ny+1:end)=eye(AuxiliaryVariables_nbr);
A11=zeros(m);A11(1:ny,1:ny)=A1;A11(1:ny,ny+1:end)=A6(:,aux);
A22=zeros(m);A22(1:ny,1:ny)=A2;
A33=zeros(m,instr_nbr);A33(1:ny,1:end)=A3;A33(ny+1:end,aux)=eye(AuxiliaryVariables_nbr);
A44=zeros(m,instr_nbr);A44(1:ny,1:end)=A4;
A55=zeros(m,exo_nbr);A55(1:ny,1:end)=A5;
WW=zeros(m);WW(1:ny,1:ny)=W;
endo_augm_id=setdiff(1:endo_nbr+AuxiliaryVariables_nbr,instr_id);

end

function v= SylvesterDoubling (d,g,h,tol,maxit)

% DOUBLES  Solves a Sylvester equation using doubling
%
%  [v,info] = doubles (g,d,h,tol,maxit)  uses a doubling algorithm
%   to solve the  Sylvester equation v  = d + g v h

v = d;
for i =1:maxit
    vadd = g*v*h;
    v = v+vadd;
    if norm (vadd,1) <= (tol*norm(v,1))
        break
    end
    g = g*g;
    h = h*h;
end

end

function v = SylvesterHessenbergSchur(d,g,h)
%
% DSYLHS  Solves a discrete time sylvester equation     using the
% Hessenberg-Schur algorithm
%
% v = DSYLHS(g,d,h) computes the matrix v that satisfies the
% discrete time sylvester equation
%
%           v = d + g'vh

if size(g,1) >= size(h,1)
    [u,gbarp] = hess(g');
    [t,hbar] = schur(h);
    [vbar] = sylvest_private(gbarp,u'*d*t,hbar,1e-15);
    v = u*vbar*t';
else
    [u,gbar] = schur(g);
    [t,hbarp] = hess(h');
    [vbar] = sylvest_private(hbarp,t'*d'*u,gbar,1e-15);
    v = u*vbar'*t';
end

end


function v = sylvest_private(g,d,h,tol)
%
% SYLVEST  Solves a Sylvester equation
%
%  solves the Sylvester equation
%         v = d + g v h
%  for v where both g and h must be  upper block triangular.
%  The output info is zero on a successful return.
%  The input tol indicates when an element of g or h should be considered
%  zero.

[m,n] = size(d);
v = zeros(m,n);
w = eye(m);
i = 1;
temp = [];

%First handle the i = 1 case outside the loop

if i< n
    if abs(h(i+1,i)) < tol
        v(:,i)= (w - g*h(i,i))\d(:,i);
        i = i+1;
    else
        A = [w-g*h(i,i)     (-g*h(i+1,i));...
             -g*h(i,i+1)    w-g*h(i+1,i+1)];
        C = [d(:,i); d(:,i+1)];
        X = A\C;
        v(:,i) = X(1:m,:);
        v(:,i+1) = X(m+1:2*m,  :);
        i = i+2;
    end
end

%Handle the rest of the matrix with the possible exception of i=n

while i<n
    b= i-1;
    temp = [temp g*v(:,size(temp,2)+1:b)]; %#ok<AGROW>
    if abs(h(i+1,i)) < tol
        v(:,i) = (w - g*h(i,i))\(d(:,i) + temp*h(1:b,i));
        i = i+1;
    else
        A = [w - g*h(i,i)    (-g*h(i+1,i));   ...
             -g*h(i,i+1)    w - g*h(i+1,i+1)];
        C = [d(:,i) + temp*h(1:b,i);         ...
             d(:,i+1) + temp*h(1:b,i+1)];
        X = A\C;
        v(:,i) = X(1:m,:);
        v(:,i+1) = X(m+1:2*m, :);
        i = i+2;
    end
end

%Handle the i = n case if i=n was not in a 2-2 block

if i==n
    b = i-1;
    temp = [temp g*v(:,size(temp,2)+1:b)];
    v(:,i) = (w-g*h(i,i))\(d(:,i) + temp*h(1:b,i));
end

end
