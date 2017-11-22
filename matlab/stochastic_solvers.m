function [dr,info] = stochastic_solvers(dr,task,M_,options_,oo_)
% function [dr,info,M_,options_,oo_] = stochastic_solvers(dr,task,M_,options_,oo_)
% computes the reduced form solution of a rational expectations model (first, second or third
% order approximation of the stochastic model around the deterministic steady state).
%
% INPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   task       [integer]          if task = 0 then dr1 computes decision rules.
%                                 if task = 1 then dr1 computes eigenvalues.
%   M_         [matlab structure] Definition of the model.
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results
%
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   info       [integer]          info=1: the model doesn't define current variables uniquely
%                                 info=2: problem in mjdgges.dll info(2) contains error code.
%                                 info=3: BK order condition not satisfied info(2) contains "distance"
%                                         absence of stable trajectory.
%                                 info=4: BK order condition not satisfied info(2) contains "distance"
%                                         indeterminacy.
%                                 info=5: BK rank condition not satisfied.
%                                 info=6: The jacobian matrix evaluated at the steady state is complex.
%                                 info=9: k_order_pert was unable to compute the solution
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   none.
%

% Copyright (C) 1996-2017 Dynare Team
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

info = 0;

if options_.linear
    options_.order = 1;
end

local_order = options_.order;
if M_.hessian_eq_zero && local_order~=1
    local_order = 1;
    warning('stochastic_solvers: using order = 1 because Hessian is equal to zero');
end
if options_.order>2 && ~options_.k_order_solver
    error('You need to set k_order_solver for order>2')
end

if (options_.aim_solver == 1) && (local_order > 1)
    error('Option "aim_solver" is incompatible with order >= 2')
end

if M_.maximum_endo_lag == 0
    if local_order >= 2
        fprintf('\nSTOCHASTIC_SOLVER: Dynare does not solve purely forward models at higher order.\n')
        fprintf('STOCHASTIC_SOLVER: To circumvent this restriction, you can add a backward-looking dummy equation of the form:\n')
        fprintf('STOCHASTIC_SOLVER: junk=0.9*junk(-1);\n')
        error(['2nd and 3rd order approximation not implemented for purely ' ...
               'forward models'])
    end
    if M_.exo_det_nbr~=0
        fprintf('\nSTOCHASTIC_SOLVER: Dynare does not solve purely forward models with var_exo_det.\n')
        fprintf('STOCHASTIC_SOLVER: To circumvent this restriction, you can add a backward-looking dummy equation of the form:\n')
        fprintf('STOCHASTIC_SOLVER: junk=0.9*junk(-1);\n')
        error(['var_exo_det not implemented for purely forward models'])
    end
end

if M_.maximum_endo_lead==0 && M_.exo_det_nbr~=0
    fprintf('\nSTOCHASTIC_SOLVER: Dynare does not solve purely backward models with var_exo_det.\n')
    fprintf('STOCHASTIC_SOLVER: To circumvent this restriction, you can add a foward-looking dummy equation of the form:\n')
    fprintf('STOCHASTIC_SOLVER: junk=0.9*junk(+1);\n')
    error(['var_exo_det not implemented for purely backwards models'])
end

if options_.k_order_solver
    if options_.risky_steadystate
        [dr,info] = dyn_risky_steadystate_solver(oo_.steady_state,M_,dr, ...
                                                 options_,oo_);
    else
        orig_order = options_.order;
        options_.order = local_order;
        dr = set_state_space(dr,M_,options_);
        [dr,info] = k_order_pert(dr,M_,options_);
        options_.order = orig_order;
    end
    return
end

klen = M_.maximum_lag + M_.maximum_lead + 1;
exo_simul = [repmat(oo_.exo_steady_state',klen,1) repmat(oo_.exo_det_steady_state',klen,1)];
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

it_ = M_.maximum_lag + 1;
z = repmat(dr.ys,1,klen);
if local_order == 1
    if (options_.bytecode)
        [chck, junk, loc_dr] = bytecode('dynamic','evaluate', z,exo_simul, ...
                                        M_.params, dr.ys, 1);
        jacobia_ = [loc_dr.g1 loc_dr.g1_x loc_dr.g1_xd];
    else
        [junk,jacobia_] = feval([M_.fname '_dynamic'],z(iyr0),exo_simul, ...
                                M_.params, dr.ys, it_);
    end
elseif local_order == 2
    if (options_.bytecode)
        [chck, junk, loc_dr] = bytecode('dynamic','evaluate', z,exo_simul, ...
                                        M_.params, dr.ys, 1);
        jacobia_ = [loc_dr.g1 loc_dr.g1_x];
    else
        [junk,jacobia_,hessian1] = feval([M_.fname '_dynamic'],z(iyr0),...
                                         exo_simul, ...
                                         M_.params, dr.ys, it_);
    end
    if options_.use_dll
        % In USE_DLL mode, the hessian is in the 3-column sparse representation
        hessian1 = sparse(hessian1(:,1), hessian1(:,2), hessian1(:,3), ...
                          size(jacobia_, 1), size(jacobia_, 2)*size(jacobia_, 2));
    end
    [infrow,infcol]=find(isinf(hessian1));
    if options_.debug
        if ~isempty(infrow)
            fprintf('\nSTOCHASTIC_SOLVER: The Hessian of the dynamic model contains Inf.\n')
            fprintf('STOCHASTIC_SOLVER: Try running model_diagnostics to find the source of the problem.\n')
            save([M_.fname '_debug.mat'],'hessian1')
        end
    end
    if ~isempty(infrow)
        info(1)=11;
        return
    end
    [nanrow,nancol]=find(isnan(hessian1));
    if options_.debug
        if ~isempty(nanrow)
            fprintf('\nSTOCHASTIC_SOLVER: The Hessian of the dynamic model contains NaN.\n')
            fprintf('STOCHASTIC_SOLVER: Try running model_diagnostics to find the source of the problem.\n')
            save([M_.fname '_debug.mat'],'hessian1')
        end
    end
    if ~isempty(nanrow)
        info(1)=12;
        return
    end
end

[infrow,infcol]=find(isinf(jacobia_));

if options_.debug
    if ~isempty(infrow)
        fprintf('\nSTOCHASTIC_SOLVER: The Jacobian of the dynamic model contains Inf. The problem is associated with:\n\n')
        display_problematic_vars_Jacobian(infrow,infcol,M_,dr.ys,'dynamic','STOCHASTIC_SOLVER: ')
        save([M_.fname '_debug.mat'],'jacobia_')
    end
end

if ~isempty(infrow)
    info(1)=10;
    return
end

if ~isreal(jacobia_)
    if max(max(abs(imag(jacobia_)))) < 1e-15
        jacobia_ = real(jacobia_);
    else
        if options_.debug
            [imagrow,imagcol]=find(abs(imag(jacobia_))>1e-15);
            fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the dynamic model contains imaginary parts. The problem arises from: \n\n')
            display_problematic_vars_Jacobian(imagrow,imagcol,M_,dr.ys,'dynamic','STOCHASTIC_SOLVER: ')
        end
        info(1) = 6;
        info(2) = sum(sum(imag(jacobia_).^2));
        return
    end
end

[nanrow,nancol]=find(isnan(jacobia_));
if options_.debug
    if ~isempty(nanrow)
        fprintf('\nSTOCHASTIC_SOLVER: The Jacobian of the dynamic model contains NaN. The problem is associated with:\n\n')
        display_problematic_vars_Jacobian(nanrow,nancol,M_,dr.ys,'dynamic','STOCHASTIC_SOLVER: ')
        save([M_.fname '_debug.mat'],'jacobia_')
    end
end

if ~isempty(nanrow)
    info(1) = 8;
    NaN_params=find(isnan(M_.params));
    info(2:length(NaN_params)+1) =  NaN_params;
    return
end

kstate = dr.kstate;
nstatic = M_.nstatic;
nfwrd = M_.nfwrd;
nspred = M_.nspred;
nboth = M_.nboth;
nsfwrd = M_.nsfwrd;
order_var = dr.order_var;
nd = size(kstate,1);
nz = nnz(M_.lead_lag_incidence);

sdyn = M_.endo_nbr - nstatic;

[junk,cols_b,cols_j] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+1, ...
                                                  order_var));
b = zeros(M_.endo_nbr,M_.endo_nbr);
b(:,cols_b) = jacobia_(:,cols_j);

if M_.maximum_endo_lead == 0
    % backward models: simplified code exist only at order == 1
    if local_order == 1
        [k1,junk,k2] = find(kstate(:,4));
        dr.ghx(:,k1) = -b\jacobia_(:,k2);
        if M_.exo_nbr
            dr.ghu =  -b\jacobia_(:,nz+1:end);
        end
        dr.eigval = eig(kalman_transition_matrix(dr,nstatic+(1:nspred),1:nspred,M_.exo_nbr));
        dr.full_rank = 1;
        if any(abs(dr.eigval) > options_.qz_criterium)
            temp = sort(abs(dr.eigval));
            nba = nnz(abs(dr.eigval) > options_.qz_criterium);
            temp = temp(nd-nba+1:nd)-1-options_.qz_criterium;
            info(1) = 3;
            info(2) = temp'*temp;
        end
    else
        fprintf('\nSTOCHASTIC_SOLVER: Dynare does not solve purely backward models at higher order.\n')
        fprintf('STOCHASTIC_SOLVER: To circumvent this restriction, you can add a forward-looking dummy equation of the form:\n')
        fprintf('STOCHASTIC_SOLVER: junk=0.9*junk(+1);\n')
        error(['2nd and 3rd order approximation not implemented for purely ' ...
               'backward models'])
    end
elseif options_.risky_steadystate
    orig_order = options_.order;
    options_.order = local_order;
    [dr,info] = dyn_risky_steadystate_solver(oo_.steady_state,M_,dr, ...
                                             options_,oo_);
    options_.order = orig_order;
else
    % If required, use AIM solver if not check only
    if (options_.aim_solver == 1) && (task == 0)
        [dr,info] = AIM_first_order_solver(jacobia_,M_,dr,options_.qz_criterium);

    else  % use original Dynare solver
        [dr,info] = dyn_first_order_solver(jacobia_,M_,dr,options_,task);
        if info(1) || task
            return
        end
    end

    if local_order > 1
        % Second order
        dr = dyn_second_order_solver(jacobia_,hessian1,dr,M_,...
                                     options_.threads.kronecker.A_times_B_kronecker_C,...
                                     options_.threads.kronecker.sparse_hessian_times_B_kronecker_C);

        % reordering second order derivatives, used for deterministic
        % variables below
        k1 = nonzeros(M_.lead_lag_incidence(:,order_var)');
        kk = [k1; length(k1)+(1:M_.exo_nbr+M_.exo_det_nbr)'];
        nk = size(kk,1);
        kk1 = reshape([1:nk^2],nk,nk);
        kk1 = kk1(kk,kk);
        hessian1 = hessian1(:,kk1(:));
    end
end


%exogenous deterministic variables
if M_.exo_det_nbr > 0
    gx = dr.gx;
    f1 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+2:end,order_var))));
    f0 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var))));
    fudet = sparse(jacobia_(:,nz+M_.exo_nbr+1:end));
    M1 = inv(f0+[zeros(M_.endo_nbr,nstatic) f1*gx zeros(M_.endo_nbr,nsfwrd-nboth)]);
    M2 = M1*f1;
    dr.ghud = cell(M_.exo_det_length,1);
    dr.ghud{1} = -M1*fudet;
    for i = 2:M_.exo_det_length
        dr.ghud{i} = -M2*dr.ghud{i-1}(end-nsfwrd+1:end,:);
    end

    if local_order > 1
        lead_lag_incidence = M_.lead_lag_incidence;
        k0 = find(lead_lag_incidence(M_.maximum_endo_lag+1,order_var)');
        k1 = find(lead_lag_incidence(M_.maximum_endo_lag+2,order_var)');
        hu = dr.ghu(nstatic+[1:nspred],:);
        hud = dr.ghud{1}(nstatic+1:nstatic+nspred,:);
        zx = [eye(nspred);dr.ghx(k0,:);gx*dr.Gy;zeros(M_.exo_nbr+M_.exo_det_nbr, ...
                                                      nspred)];
        zu = [zeros(nspred,M_.exo_nbr); dr.ghu(k0,:); gx*hu; zeros(M_.exo_nbr+M_.exo_det_nbr, ...
                                                          M_.exo_nbr)];
        zud=[zeros(nspred,M_.exo_det_nbr);dr.ghud{1};gx(:,1:nspred)*hud;zeros(M_.exo_nbr,M_.exo_det_nbr);eye(M_.exo_det_nbr)];
        R1 = hessian1*kron(zx,zud);
        dr.ghxud = cell(M_.exo_det_length,1);
        kf = [M_.endo_nbr-nfwrd-nboth+1:M_.endo_nbr];
        kp = nstatic+[1:nspred];
        dr.ghxud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{1}(kp,:)));
        Eud = eye(M_.exo_det_nbr);
        for i = 2:M_.exo_det_length
            hudi = dr.ghud{i}(kp,:);
            zudi=[zeros(nspred,M_.exo_det_nbr);dr.ghud{i};gx(:,1:nspred)*hudi;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
            R2 = hessian1*kron(zx,zudi);
            dr.ghxud{i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(dr.Gy,Eud)+dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{i}(kp,:)))-M1*R2;
        end
        R1 = hessian1*kron(zu,zud);
        dr.ghudud = cell(M_.exo_det_length,1);
        dr.ghuud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghu(kp,:),dr.ghud{1}(kp,:)));
        Eud = eye(M_.exo_det_nbr);
        for i = 2:M_.exo_det_length
            hudi = dr.ghud{i}(kp,:);
            zudi=[zeros(nspred,M_.exo_det_nbr);dr.ghud{i};gx(:,1:nspred)*hudi;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
            R2 = hessian1*kron(zu,zudi);
            dr.ghuud{i} = -M2*dr.ghxud{i-1}(kf,:)*kron(hu,Eud)-M1*R2;
        end
        R1 = hessian1*kron(zud,zud);
        dr.ghudud = cell(M_.exo_det_length,M_.exo_det_length);
        dr.ghudud{1,1} = -M1*R1-M2*dr.ghxx(kf,:)*kron(hud,hud);
        for i = 2:M_.exo_det_length
            hudi = dr.ghud{i}(nstatic+1:nstatic+nspred,:);
            zudi=[zeros(nspred,M_.exo_det_nbr);dr.ghud{i};gx(:,1:nspred)*hudi+dr.ghud{i-1}(kf,:);zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
            R2 = hessian1*kron(zudi,zudi);
            dr.ghudud{i,i} = -M2*(dr.ghudud{i-1,i-1}(kf,:)+...
                                  2*dr.ghxud{i-1}(kf,:)*kron(hudi,Eud) ...
                                  +dr.ghxx(kf,:)*kron(hudi,hudi))-M1*R2;
            R2 = hessian1*kron(zud,zudi);
            dr.ghudud{1,i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hud,Eud)+...
                                  dr.ghxx(kf,:)*kron(hud,hudi))...
                -M1*R2;
            for j=2:i-1
                hudj = dr.ghud{j}(kp,:);
                zudj=[zeros(nspred,M_.exo_det_nbr);dr.ghud{j};gx(:,1:nspred)*hudj;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
                R2 = hessian1*kron(zudj,zudi);
                dr.ghudud{j,i} = -M2*(dr.ghudud{j-1,i-1}(kf,:)+dr.ghxud{j-1}(kf,:)* ...
                                      kron(hudi,Eud)+dr.ghxud{i-1}(kf,:)* ...
                                      kron(hudj,Eud)+dr.ghxx(kf,:)*kron(hudj,hudi))-M1*R2;
            end
        end
    end
end

if options_.loglinear
    % this needs to be extended for order=2,3
    [il,il1,ik,k1] = indices_lagged_leaded_exogenous_variables(dr.order_var,M_);
    [illag,illag1,iklag,klag1] = indices_lagged_leaded_exogenous_variables(dr.order_var(M_.nstatic+(1:M_.nspred)),M_);
    if ~isempty(ik)
        if M_.nspred > 0
            dr.ghx(ik,iklag) = repmat(1./dr.ys(k1),1,length(klag1)).*dr.ghx(ik,iklag).* ...
                repmat(dr.ys(klag1)',length(ik),1);
            dr.ghx(ik,illag) = repmat(1./dr.ys(k1),1,length(illag)).*dr.ghx(ik,illag);
        end
        if M_.exo_nbr > 0
            dr.ghu(ik,:) = repmat(1./dr.ys(k1),1,M_.exo_nbr).*dr.ghu(ik,:);
        end
    end
    if ~isempty(il) && M_.nspred > 0
        dr.ghx(il,iklag) = dr.ghx(il,iklag).*repmat(dr.ys(klag1)', ...
                                                    length(il),1);
    end
    if local_order > 1
        error('Loglinear options currently only works at order 1')
    end
end
end
