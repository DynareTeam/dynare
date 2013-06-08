function y = solve_two_boundaries(fname, y, x, params, steady_state, y_index, nze, periods, y_kmin_l, y_kmax_l, is_linear, Block_Num, y_kmin, maxit_, solve_tolf, lambda, cutoff, stack_solve_algo, M, oo)
% Computes the deterministic simulation of a block of equation containing
% both lead and lag variables using relaxation methods 
%
% INPUTS
%   fname               [string]        name of the file containing the block
%                                       to simulate
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   params              [vector]        All the parameters of the model
%   steady_state        [vector]        steady state of the model
%   y_index             [vector of int] The index of the endogenous variables of
%                                       the block
%   nze                 [integer]       number of non-zero elements in the
%                                       jacobian matrix
%   periods             [integer]       number of simulation periods
%   y_kmin_l            [integer]       maximum number of lag in the block
%   y_kmax_l            [integer]       maximum number of lead in the block
%   is_linear           [integer]       if is_linear=1 the block is linear
%                                       if is_linear=0 the block is not linear
%   Block_Num           [integer]       block number
%   y_kmin              [integer]       maximum number of lag in the model
%   maxit_              [integer]       maximum number of iteration in Newton
%   solve_tolf          [double]        convergence criteria
%   lambda              [double]        initial value of step size in
%   Newton
%   cutoff              [double]        cutoff to correct the direction in Newton in case
%                                       of singular jacobian matrix
%   stack_solve_algo    [integer]       linear solver method used in the
%                                       Newton algorithm : 
%                                            - 1 sprse LU
%                                            - 2 GMRES
%                                            - 3 BicGStab
%                                            - 4 Optimal path length
%
% OUTPUTS
%   y                   [matrix]        All endogenous variables of the model      
%  
% ALGORITHM
%   Newton with LU or GMRES or BicGstab
%    
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 1996-2013 Dynare Team
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

cvg=0;
iter=0;
Per_u_=0;
g2 = [];
g3 = [];
Blck_size=size(y_index,2);
correcting_factor=0.01;
ilu_setup.droptol=1e-10;
ilu_setup.type = 'ilutp';
%ilu_setup.milu = 'col';
ilu_setup.milu = 'off';
ilu_setup.thresh = 1;
ilu_setup.udiag = 0;
max_resa=1e100;
Jacobian_Size=Blck_size*(y_kmin+y_kmax_l +periods);
g1=spalloc( Blck_size*periods, Jacobian_Size, nze*periods);
reduced = 0;
while ~(cvg==1 || iter>maxit_),
    [r, y, g1, g2, g3, b]=feval(fname, y, x, params, steady_state, periods, 0, y_kmin, Blck_size);
    preconditioner = 2;
    g1a=g1(:, y_kmin*Blck_size+1:(periods+y_kmin)*Blck_size);
    term1 = g1(:, 1:y_kmin_l*Blck_size)*reshape(y(1+y_kmin-y_kmin_l:y_kmin,y_index)',1,y_kmin_l*Blck_size)';
    term2 = g1(:, (periods+y_kmin_l)*Blck_size+1:(periods+y_kmin_l+y_kmax_l)*Blck_size)*reshape(y(periods+y_kmin+1:periods+y_kmin+y_kmax_l,y_index)',1,y_kmax_l*Blck_size)';
    b = b - term1 - term2;
    [max_res, max_indx]=max(max(abs(r')));
    if(~isreal(r))
        max_res = (-max_res^2)^0.5;
    end;
    %     if(~isreal(r))
    %       max_res=(-(max(max(abs(r))))^2)^0.5;
    %     else
    %       max_res=max(max(abs(r)));
    %     end;
    if(~isreal(max_res) || isnan(max_res))
        cvg = 0;
    elseif(is_linear && iter>0)
        cvg = 1;
    else
        cvg=(max_res<solve_tolf);
    end;
    if(~cvg)
        if(iter>0)
            if(~isreal(max_res) || isnan(max_res) || (max_resa<max_res && iter>1))
                if(~isreal(max_res))
                    disp(['Variable ' M.endo_names(max_indx,:) ' (' int2str(max_indx) ') returns an undefined value']);
                end;
                if(isnan(max_res))
                    detJ=det(g1aa);
                    if(abs(detJ)<1e-7)
                        max_factor=max(max(abs(g1aa)));
                        ze_elem=sum(diag(g1aa)<cutoff);
                        disp([num2str(full(ze_elem),'%d') ' elements on the Jacobian diagonal are below the cutoff (' num2str(cutoff,'%f') ')']);
                        if(correcting_factor<max_factor)
                            correcting_factor=correcting_factor*4;
                            disp(['The Jacobain matrix is singular, det(Jacobian)=' num2str(detJ,'%f') '.']);
                            disp(['    trying to correct the Jacobian matrix:']);
                            disp(['    correcting_factor=' num2str(correcting_factor,'%f') ' max(Jacobian)=' num2str(full(max_factor),'%f')]);
                            dx = (g1aa+correcting_factor*speye(periods*Blck_size))\ba- ya;
                            y(1+y_kmin:periods+y_kmin,y_index)=reshape((ya_save+lambda*dx)',length(y_index),periods)';
                            continue;
                        else
                            disp('The singularity of the jacobian matrix could not be corrected');
                            return;
                        end;
                    end;
                elseif(lambda>1e-8)
                    lambda=lambda/2;
                    reduced = 1;
                    disp(['reducing the path length: lambda=' num2str(lambda,'%f')]);
                    y(1+y_kmin:periods+y_kmin,y_index)=reshape((ya_save+lambda*dx)',length(y_index),periods)';
                    continue;
                else
                    if(cutoff == 0)
                        fprintf('Error in simul: Convergence not achieved in block %d, after %d iterations.\n Increase "options_.maxit_".\n',Block_Num, iter);
                    else
                        fprintf('Error in simul: Convergence not achieved in block %d, after %d iterations.\n Increase "options_.maxit_" or set "cutoff=0" in model options.\n',Block_Num, iter);
                    end;
                    oo.deterministic_simulation.status = 0;
                    oo.deterministic_simulation.error = max_res;
                    oo.deterministic_simulation.iterations = iter;
                    oo.deterministic_simulation.block(Block_Num).status = 0;% Convergency failed.
                    oo.deterministic_simulation.block(Block_Num).error = max_res;
                    oo.deterministic_simulation.block(Block_Num).iterations = iter;
                    return;
                end;
            else
                if(lambda<1)
                    lambda=max(lambda*2, 1);
                end;
            end;
        end;
        ya = reshape(y(y_kmin+1:y_kmin+periods,y_index)',1,periods*Blck_size)';
        ya_save=ya;
        g1aa=g1a;
        ba=b;
        max_resa=max_res;
        if(stack_solve_algo==0),
            dx = g1a\b- ya;
            ya = ya + lambda*dx;
            y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';
        elseif(stack_solve_algo==1),
            for t=1:periods;
                first_elem = (t-1)*Blck_size+1;
                last_elem = t*Blck_size;
                next_elem = (t+1)*Blck_size;
                Elem = first_elem:last_elem;
                Elem_1 = last_elem+1:next_elem;
                B1_inv = inv(g1a(Elem, Elem));
                if (t < periods)
                    S1 = B1_inv * g1a(Elem, Elem_1);
                end;
                g1a(Elem, Elem_1) = S1;
                b(Elem) = B1_inv * b(Elem);
                g1a(Elem, Elem) = ones(Blck_size, Blck_size);
                if (t < periods)
                    g1a(Elem_1, Elem_1) = g1a(Elem_1, Elem_1) - g1a(Elem_1, Elem) * S1;
                    b(Elem_1) = b(Elem_1) - g1a(Elem_1, Elem) * b(Elem);
                    g1a(Elem_1, Elem) = zeros(Blck_size, Blck_size);
                end;
            end;
            za = b(Elem);
            zaa = za;
            y_Elem = (periods - 1) * Blck_size + 1:(periods) * Blck_size;
            dx = ya;
            dx(y_Elem) = za - ya(y_Elem);
            ya(y_Elem) = ya(y_Elem) + lambda*dx(y_Elem);
            for t=periods-1:-1:1;
                first_elem = (t-1)*Blck_size+1;
                last_elem = t*Blck_size;
                next_elem = (t+1)*Blck_size;
                Elem_1 = last_elem+1:next_elem;
                Elem = first_elem:last_elem;
                za = b(Elem) - g1a(Elem, Elem_1) * zaa;
                zaa = za;
                y_Elem = Blck_size * (t-1)+1:Blck_size * (t);
                dx(y_Elem) = za - ya(y_Elem);
                ya(y_Elem) = ya(y_Elem) + lambda*dx(y_Elem);
                y(y_kmin + t, y_index) = ya(y_Elem);
            end;
        elseif(stack_solve_algo==2),
            flag1=1;
            while(flag1>0)
                if preconditioner == 2
                    [L1, U1]=ilu(g1a,ilu_setup);
                elseif preconditioner == 3
                    Size = Blck_size;
                     gss1 =  g1a(Size + 1: 2*Size,Size + 1: 2*Size) + g1a(Size + 1: 2*Size,2*Size+1: 3*Size);
                     [L1, U1]=lu(gss1);
                     L(1:Size,1:Size) = L1;
                     U(1:Size,1:Size) = U1;
                     gss2 = g1a(Size + 1: 2*Size,1: Size) + g1a(Size + 1: 2*Size,Size+1: 2*Size) + g1a(Size + 1: 2*Size,2*Size+1: 3*Size);
                     [L2, U2]=lu(gss2);
                     L(Size+1:(periods-1)*Size,Size+1:(periods-1)*Size) = kron(eye(periods-2), L2);
                     U(Size+1:(periods-1)*Size,Size+1:(periods-1)*Size) = kron(eye(periods-2), U2);
                     gss2 = g1a(Size + 1: 2*Size,1: Size) + g1a(Size + 1: 2*Size,Size+1: 2*Size);
                     [L3, U3]=lu(gss2);
                     L((periods-1)*Size+1:periods*Size,(periods-1)*Size+1:periods*Size) = L3;
                     U((periods-1)*Size+1:periods*Size,(periods-1)*Size+1:periods*Size) = U3;
                     L1 = L;
                     U1 = U;
                elseif preconditioner == 4
                    Size = Blck_size;
                    gss1 =  g1a(1: 3*Size, 1: 3*Size);
                    [L, U] = lu(gss1);
                    L1 = kron(eye(ceil(periods/3)),L);
                    U1 = kron(eye(ceil(periods/3)),U);
                    L1 = L1(1:periods * Size, 1:periods * Size);
                    U1 = U1(1:periods * Size, 1:periods * Size);
                end;
                [za,flag1] = gmres(g1a,b,Blck_size,1e-6,Blck_size*periods,L1,U1);
                if (flag1>0 || reduced)
                    if(flag1==1)
                        disp(['Error in simul: No convergence inside GMRES after ' num2str(periods*10,'%6d') ' iterations, in block ' num2str(Blck_size,'%3d')]);
                    elseif(flag1==2)
                        disp(['Error in simul: Preconditioner is ill-conditioned, in block ' num2str(Blck_size,'%3d')]);
                    elseif(flag1==3)
                        disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block ' num2str(Blck_size,'%3d')]);
                    end;
                    ilu_setup.droptol = ilu_setup.droptol/10;
                    reduced = 0;
                else
                    dx = za - ya;
                    ya = ya + lambda*dx;
                    y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';
                end;
            end;
        elseif(stack_solve_algo==3),
            flag1=1;
            while(flag1>0)
                if preconditioner == 2
                    [L1, U1]=ilu(g1a,ilu_setup);
                    [za,flag1] = bicgstab(g1a,b,1e-7,Blck_size*periods,L1,U1);
                elseif (preconditioner == 3)
                    Size = Blck_size;
                    gss0 = g1a(Size + 1: 2*Size,1: Size) + g1a(Size + 1: 2*Size,Size+1: 2*Size) + g1a(Size + 1: 2*Size,2*Size+1: 3*Size);
                    [L1, U1]=lu(gss0);
                    P1 = eye(size(gss0));
                    Q1 = eye(size(gss0));
                    L = kron(eye(periods),L1);
                    U = kron(eye(periods),U1);
                    P = kron(eye(periods),P1);
                    Q = kron(eye(periods),Q1);
                    [za,flag1] = bicgstab1(g1a,b,1e-7,Blck_size*periods,L,U, P, Q);
                else
                    Size = Blck_size;
                    gss0 = g1a(Size + 1: 2*Size,1: Size) + g1a(Size + 1: 2*Size,Size+1: 2*Size) + g1a(Size + 1: 2*Size,2*Size+1: 3*Size);
                    [L1, U1]=lu(gss0);
                    L1 = kron(eye(periods),L1);
                    U1 = kron(eye(periods),U1);
                    [za,flag1] = bicgstab(g1a,b,1e-7,Blck_size*periods,L1,U1);
                end;
                if (flag1>0 || reduced)
                    if(flag1==1)
                        disp(['Error in simul: No convergence inside BICGSTAB after ' num2str(periods*10,'%6d') ' iterations, in block ' num2str(Blck_size,'%3d')]);
                    elseif(flag1==2)
                        disp(['Error in simul: Preconditioner is ill-conditioned, in block ' num2str(Blck_size,'%3d')]);
                    elseif(flag1==3)
                        disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block ' num2str(Blck_size,'%3d')]);
                    end;
                    ilu_setup.droptol = ilu_setup.droptol/10;
                    reduced = 0;
                else
                    dx = za - ya;
                    ya = ya + lambda*dx;
                    y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';
                end;
            end;
        elseif(stack_solve_algo==4),
            ra = reshape(r(:, y_kmin+1:periods+y_kmin),periods*Blck_size, 1);
            stpmx = 100 ;
            stpmax = stpmx*max([sqrt(ya'*ya);size(y_index,2)]);
            nn=1:size(ra,1);
            g = (ra'*g1a)';
            f = 0.5*ra'*ra;
            p = -g1a\ra;
            [yn,f,ra,check]=lnsrch1(ya,f,g,p,stpmax, ...
                                    'lnsrch1_wrapper_two_boundaries',nn,nn,  fname, y, y_index, x, params, steady_state, periods, y_kmin, Blck_size);
            dx = ya - yn;
            y(1+y_kmin:periods+y_kmin,y_index)=reshape(yn',length(y_index),periods)';
        end
    end
    iter=iter+1;
    disp(['iteration: ' num2str(iter,'%d') ' error: ' num2str(max_res,'%e')]);
end;
if (iter>maxit_)
    disp(['No convergence after ' num2str(iter,'%4d') ' iterations in Block ' num2str(Block_Num,'%d')]);
    oo.deterministic_simulation.status = 0;
    oo.deterministic_simulation.error = max_res;
    oo.deterministic_simulation.iterations = iter;
    oo.deterministic_simulation.block(Block_Num).status = 0;% Convergency failed.
    oo.deterministic_simulation.block(Block_Num).error = max_res;
    oo.deterministic_simulation.block(Block_Num).iterations = iter;
    return;
end
oo.deterministic_simulation.status = 1;
oo.deterministic_simulation.error = max_res;
oo.deterministic_simulation.iterations = iter;
oo.deterministic_simulation.block(Block_Num).status = 1;% Convergency obtained.
oo.deterministic_simulation.block(Block_Num).error = max_res;
oo.deterministic_simulation.block(Block_Num).iterations = iter;
return;
