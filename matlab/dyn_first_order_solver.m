function [dr,info] = dyn_first_order_solver(jacobia,DynareModel,dr,DynareOptions,task)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info}] =} dyn_first_order_solver (@var{jacobia},@var{DynareModel},@var{dr},@var{DynareOptions},@var{task})
%! @anchor{dyn_first_order_solver}
%! @sp 1
%! Computes the first order reduced form of the DSGE model
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item jacobia
%! Matrix containing the Jacobian of the model
%! @item DynareModel
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item qz_criterium
%! Double containing the criterium to separate explosive from stable eigenvalues
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==7
%! One of the generalized eigenvalues is close to 0/0
%! @end table
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2012 Dynare Team
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

persistent reorder_jacobian_columns innovations_idx index_s index_m index_c
persistent index_p row_indx index_0m index_0p k1 k2 state_var
persistent ndynamic nstatic nfwrd npred nboth nd nsfwrd n_current index_d 
persistent index_e index_d1 index_d2 index_e1 index_e2 row_indx_de_1 
persistent row_indx_de_2 cols_b


if ~nargin
    if nargout
        error('dyn_first_order_solver:: Initialization mode returns zero argument!')
    end
    reorder_jacobian_columns = [];
    return
end

exo_nbr = DynareModel.exo_nbr;

if isempty(reorder_jacobian_columns)
    
    maximum_lag = DynareModel.maximum_endo_lag;
    kstate   = dr.kstate;
    nfwrd    = DynareModel.nfwrd;
    nboth    = DynareModel.nboth;
    npred    = DynareModel.npred;
    nstatic  = DynareModel.nstatic;
    ndynamic = DynareModel.ndynamic;
    nsfwrd   = DynareModel.nsfwrd;
    n        = DynareModel.endo_nbr;

    k1 = 1:(npred+nboth);
    k2 = 1:(nfwrd+nboth);

    order_var = dr.order_var;
    nd = size(kstate,1);
    lead_lag_incidence = DynareModel.lead_lag_incidence;
    nz = nnz(lead_lag_incidence);

    lead_id = find(lead_lag_incidence(maximum_lag+2,:));
    lead_idx = lead_lag_incidence(maximum_lag+2,lead_id);
    if maximum_lag
        lag_id = find(lead_lag_incidence(1,:));
        lag_idx = lead_lag_incidence(1,lag_id);
        static_id = find((lead_lag_incidence(1,:) == 0) & ...
                         (lead_lag_incidence(3,:) == 0));
    else
        lag_id = [];
        lag_idx = [];
        static_id = find(lead_lag_incidence(2,:)==0);
    end

    both_id = intersect(lead_id,lag_id);
    if maximum_lag
        no_both_lag_id = setdiff(lag_id,both_id);
    else
        no_both_lag_id = lag_id;
    end
    innovations_idx = nz+(1:exo_nbr);
    state_var  = [no_both_lag_id, both_id];

    index_c  = nonzeros(lead_lag_incidence(maximum_lag+1,:));             % Index of all endogenous variables present at time=t
    n_current = length(index_c);

    index_s  = npred+nboth+(1:nstatic);     % Index of all static
                                            % endogenous variables
                                            % present at time=t
    indexi_0 = npred+nboth;

    npred0 = nnz(lead_lag_incidence(maximum_lag+1,no_both_lag_id));
    index_0m = indexi_0+nstatic+(1:npred0);
    nfwrd0 = nnz(lead_lag_incidence(2,lead_id));
    index_0p = indexi_0+nstatic+npred0+(1:nfwrd0);
    index_m  = 1:(npred+nboth);
    index_p  = npred+nboth+n_current+(1:nfwrd+nboth);
    row_indx_de_1 = 1:ndynamic;
    row_indx_de_2 = ndynamic+(1:nboth);
    row_indx = nstatic+row_indx_de_1;
    index_d = [index_0m index_p];
    llx = lead_lag_incidence(:,order_var);
    index_d1 = [find(llx(maximum_lag+1,nstatic+(1:npred))), npred+nboth+(1:nsfwrd) ];
    index_d2 = npred+(1:nboth);

    index_e = [index_m index_0p];
    index_e1 = [1:npred+nboth, npred+nboth+find(llx(maximum_lag+1,nstatic+npred+(1: ...
                                                      nsfwrd)))];
    index_e2 = npred+nboth+(1:nboth);
    
    [junk,cols_b] = find(lead_lag_incidence(maximum_lag+1, order_var));

    reorder_jacobian_columns = [nonzeros(lead_lag_incidence(:,order_var)'); ...
                        nz+(1:exo_nbr)'];
end

dr.ghx = [];
dr.ghu = [];

dr.state_var = state_var;

jacobia = jacobia(:,reorder_jacobian_columns);

if nstatic > 0
    [Q, junk] = qr(jacobia(:,index_s));
    aa = Q'*jacobia;
else
    aa = jacobia;
end

A = aa(:,index_m);  % Jacobain matrix for lagged endogeneous variables
B(:,cols_b) = aa(:,index_c);  % Jacobian matrix for contemporaneous endogeneous variables
C = aa(:,index_p);  % Jacobain matrix for led endogeneous variables

info = 0;
if task ~= 1 && (DynareOptions.dr_cycle_reduction || DynareOptions.dr_logarithmic_reduction)
    if n_current < DynareModel.endo_nbr
        if DynareOptions.dr_cycle_reduction
            error(['The cycle reduction algorithme can''t be used when the ' ...
               'coefficient matrix for current variables isn''t invertible'])
        elseif DynareOptions.dr_logarithmic_reduction
            error(['The logarithmic reduction algorithme can''t be used when the ' ...
                   'coefficient matrix for current variables isn''t invertible'])
        end
    end  
    A1 = [aa(row_indx,index_m ) zeros(ndynamic,nfwrd)];
    B1 = [aa(row_indx,index_0m) aa(row_indx,index_0p) ];
    C1 = [zeros(ndynamic,npred) aa(row_indx,index_p)];
    if DynareOptions.gpu
        A1 = gpuArray(A1);
        B1 = gpuArray(B1);
        C1 = gpuArray(C1);
    end
    if DynareOptions.dr_cycle_reduction == 1
        [ghx, info] = cycle_reduction(A1, B1, C1, DynareOptions.dr_cycle_reduction_tol);
    else
        [ghx, info] = logarithmic_reduction(C1, B1, A1, DynareOptions.dr_logarithmic_reduction_tol, DynareOptions.dr_logarithmic_reduction_maxiter);
    end
    if info
        % cycle_reduction or logarithmic redution failed and set info
        return
    end
    ghx = ghx(:,index_m);
    hx = ghx(1:npred+nboth,:);
    gx = ghx(1+npred:end,:);
else
    D = zeros(ndynamic+nboth);
    E = zeros(ndynamic+nboth);
    D(row_indx_de_1,index_d1) = aa(row_indx,index_d);
    D(row_indx_de_2,index_d2) = eye(nboth);
    E(row_indx_de_1,index_e1) = -aa(row_indx,index_e);
    E(row_indx_de_2,index_e2) = eye(nboth);
    
    [err, ss, tt, w, sdim, dr.eigval, info1] = mjdgges(E, D, DynareOptions.qz_criterium, DynareOptions.qz_zero_threshold);
    mexErrCheck('mjdgges', err);

    if info1
        if info1 == -30
            % one eigenvalue is close to 0/0
            info(1) = 7;
        else
            info(1) = 2;
            info(2) = info1;
            info(3) = size(E,2);
        end
        return
    end

    nba = nd-sdim;

    if task == 1
        if rcond(w(npred+nboth+1:end,npred+nboth+1:end)) < 1e-9
            dr.full_rank = 0;
        else
            dr.full_rank = 1;
        end            
        return
    end

    if nba ~= nsfwrd
        temp = sort(abs(dr.eigval));
        if nba > nsfwrd
            temp = temp(nd-nba+1:nd-nsfwrd)-1-DynareOptions.qz_criterium;
            info(1) = 3;
        elseif nba < nsfwrd;
            temp = temp(nd-nsfwrd+1:nd-nba)-1-DynareOptions.qz_criterium;
            info(1) = 4;
        end
        info(2) = temp'*temp;
        return
    end

    %First order approximation
    indx_stable_root = 1: (nd - nsfwrd);     %=> index of stable roots
    indx_explosive_root = npred + nboth + 1:nd;  %=> index of explosive roots
                                                 % derivatives with respect to dynamic state variables
                                                 % forward variables
    Z = w';
    Z11 = Z(indx_stable_root,    indx_stable_root);
    Z21  = Z(indx_explosive_root, indx_stable_root);
    Z22  = Z(indx_explosive_root, indx_explosive_root);
    [minus_gx,rc] = linsolve(Z22,Z21);
    if rc < 1e-9
        % Z22 is near singular
        info(1) = 5;
        info(2) = -log(rc);
        return;
    end
    gx  = -minus_gx;
    % predetermined variables
    opts.UT = true;
    opts.TRANSA = true;
    hx1 = linsolve(tt(indx_stable_root, indx_stable_root),Z11,opts)';
    hx2 = linsolve(Z11,ss(indx_stable_root, indx_stable_root)')';
    hx =  hx1*hx2;
    ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];
end

if nstatic > 0
    B_static = B(:,1:nstatic);  % submatrix containing the derivatives w.r. to static variables
else
    B_static = [];
end;
%static variables, backward variable, mixed variables and forward variables
B_pred = B(:,nstatic+1:nstatic+npred+nboth);
B_fyd = B(:,nstatic+npred+nboth+1:end);

% static variables
if nstatic > 0
    temp = - C(1:nstatic,:)*gx*hx;
    b(:,cols_b) = aa(:,index_c);
    b10 = b(1:nstatic, 1:nstatic);
    b11 = b(1:nstatic, nstatic+1:end);
    temp(:,index_m) = temp(:,index_m)-A(1:nstatic,:);
    temp = b10\(temp-b11*ghx);
    ghx = [temp; ghx];
    temp = [];
end

A_ = real([B_static C*gx+B_pred B_fyd]); % The state_variable of the block are located at [B_pred B_both]

if exo_nbr
    if nstatic > 0
        fu = Q' * jacobia(:,innovations_idx);
    else
        fu = jacobia(:,innovations_idx);
    end;

    ghu = - A_ \ fu;
else
    ghu = [];
end;

dr.ghx = ghx;
dr.ghu = ghu;

if DynareOptions.aim_solver ~= 1 && DynareOptions.use_qzdiv
    % Necessary when using Sims' routines for QZ
    dr.ghx = real(ghx);
    dr.ghu = real(ghu);
    hx = real(hx);
end

% non-predetermined variables
dr.gx = gx;
%predetermined (endogenous state) variables, square transition matrix
dr.Gy = hx;
