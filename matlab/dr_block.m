function [dr,info,M_,options_,oo_] = dr_block(dr,task,M_,options_,oo_)
% function [dr,info,M_,options_,oo_] = dr_block(dr,task,M_,options_,oo_)
% computes the reduced form solution of a rational expectation model (first
% approximation of the stochastic model around the deterministic steady state). 
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
%   M_         [matlab structure]            
%   options_   [matlab structure]
%   oo_        [matlab structure]
%  
% ALGORITHM
%   first order block relaxation method applied to the model
%    E[A Yt-1 + B Yt + C Yt+1 + ut] = 0
%    
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 2010-2012 Dynare Team
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
verbose = 0;
if options_.order > 1
    error('2nd and 3rd order approximation not implemented with block option')
end

z = repmat(dr.ys,1,M_.maximum_lead + M_.maximum_lag + 1);
zx = repmat([oo_.exo_simul oo_.exo_det_simul],M_.maximum_lead + M_.maximum_lag + 1, 1);
if (isfield(M_,'block_structure'))
    data = M_.block_structure.block;
    Size = length(M_.block_structure.block);
else
    data = M_;
    Size = 1;
end;
if (options_.bytecode)
    [chck, zz, data]= bytecode('dynamic','evaluate', z, zx, M_.params, dr.ys, 1, data);
else
    [r, data] = feval([M_.fname '_dynamic'], z', zx, M_.params, dr.ys, M_.maximum_lag+1, data);
    chck = 0;
end;
mexErrCheck('bytecode', chck);
dr.full_rank = 1;
dr.eigval = [];
dr.nd = 0;

dr.ghx = [];
dr.ghu = [];
%Determine the global list of state variables:
dr.state_var = M_.state_var;
M_.block_structure.state_var = dr.state_var;
n_sv = size(dr.state_var, 2);
dr.ghx = zeros(M_.endo_nbr, length(dr.state_var));
dr.exo_var = 1:M_.exo_nbr;
dr.ghu = zeros(M_.endo_nbr, M_.exo_nbr);
for i = 1:Size;
    ghx = [];
    indexi_0 = 0;
    if (verbose)
        disp('======================================================================');
        disp(['Block ' int2str(i)]);
        disp('-----------');
        data(i)
    end;
    n_pred = data(i).n_backward;
    n_fwrd = data(i).n_forward;
    n_both = data(i).n_mixed;
    n_static = data(i).n_static;
    nd = n_pred + n_fwrd + 2*n_both;
    dr.nd = dr.nd + nd;
    n_dynamic = n_pred + n_fwrd + n_both;
    exo_nbr = M_.block_structure.block(i).exo_nbr;
    exo_det_nbr = M_.block_structure.block(i).exo_det_nbr;
    other_endo_nbr = M_.block_structure.block(i).other_endo_nbr;
    jacob = full(data(i).g1);
    lead_lag_incidence = data(i).lead_lag_incidence;
    endo = data(i).variable;
    exo = data(i).exogenous;
    if (verbose)
        disp('jacob');
        disp(jacob);
        disp('lead_lag_incidence');
        disp(lead_lag_incidence);
    end;
    maximum_lag = data(i).maximum_endo_lag;
    maximum_lead = data(i).maximum_endo_lead;
    n = n_dynamic + n_static;
    
    block_type = M_.block_structure.block(i).Simulation_Type;
    if task ~= 1
        if block_type == 2 || block_type == 4 || block_type == 7 
            block_type = 8;
        end;
    end;
    if maximum_lag > 0 && (n_pred > 0  || n_both > 0) && block_type ~= 1 
        indexi_0 = min(lead_lag_incidence(2,:));
    end;
    switch block_type
      case 1
      %% ------------------------------------------------------------------
        %Evaluate Forward
        if maximum_lag > 0 && n_pred > 0
            indx_r = find(M_.block_structure.block(i).lead_lag_incidence(1,:));
            indx_c = M_.block_structure.block(i).lead_lag_incidence(1,indx_r);
            data(i).eigval = diag(jacob(indx_r, indx_c));
            data(i).rank = 0;
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end
        dr.eigval = [dr.eigval ; data(i).eigval];
        %First order approximation
        if task ~= 1
            [tmp1, tmp2, indx_c] = find(M_.block_structure.block(i).lead_lag_incidence(2,:));
            B = jacob(:,indx_c);
            if (maximum_lag > 0 && n_pred > 0)
                [indx_r, tmp1, indx_r_v]  = find(M_.block_structure.block(i).lead_lag_incidence(1,:));
                ghx = - B \ jacob(:,indx_r_v);
            end;
            if other_endo_nbr
                fx = data(i).g1_o;
                % retrieves the derivatives with respect to endogenous
                % variable belonging to previous blocks
                fx_tm1 = zeros(n,other_endo_nbr);
                fx_t = zeros(n,other_endo_nbr);
                fx_tp1 = zeros(n,other_endo_nbr);
                % stores in fx_tm1 the lagged values of fx
                [r, c, lag] = find(data(i).lead_lag_incidence_other(1,:));
                fx_tm1(:,c) = fx(:,lag);
                % stores in fx the current values of fx
                [r, c, curr] = find(data(i).lead_lag_incidence_other(2,:));
                fx_t(:,c) = fx(:,curr);
                % stores in fx_tp1 the leaded values of fx
                [r, c, lead] = find(data(i).lead_lag_incidence_other(3,:));
                fx_tp1(:,c) = fx(:,lead);

                l_x = dr.ghx(data(i).other_endogenous,:);
                l_x_sv = dr.ghx(dr.state_var, 1:n_sv);

                selector_tm1 = M_.block_structure.block(i).tm1;
               
                ghx_other = - B \ (fx_t * l_x + (fx_tp1 * l_x * l_x_sv) + fx_tm1 * selector_tm1);
                dr.ghx(endo, :) = dr.ghx(endo, :) + ghx_other;
            end;
            
            if exo_nbr
                fu = data(i).g1_x;
                exo = dr.exo_var;
                if other_endo_nbr > 0
                    l_u_sv = dr.ghu(dr.state_var,:);
                    l_x = dr.ghx(data(i).other_endogenous,:);
                    l_u = dr.ghu(data(i).other_endogenous,:);
                    fu_complet = zeros(n, M_.exo_nbr);
                    fu_complet(:,data(i).exogenous) = fu;
                    ghu = - B \ (fu_complet + fx_tp1 * l_x * l_u_sv + (fx_t) * l_u );
                else
                    fu_complet = zeros(n, M_.exo_nbr);
                    fu_complet(:,data(i).exogenous) = fu;
                    ghu = - B \ fu_complet;
                end;
            else
                exo = dr.exo_var;
                if other_endo_nbr > 0
                     l_u_sv = dr.ghu(dr.state_var,:);
                     l_x = dr.ghx(data(i).other_endogenous,:);
                     l_u = dr.ghu(data(i).other_endogenous,:);
                     ghu = -B \ (fx_tp1 * l_x * l_u_sv + (fx_t) * l_u );
                else
                    ghu = [];
                end
            end
        end
      case 2
      %% ------------------------------------------------------------------
        %Evaluate Backward
        if maximum_lead > 0 && n_fwrd > 0
            indx_r = find(M_.block_structure.block(i).lead_lag_incidence(3,:));
            indx_c = M_.block_structure.block(i).lead_lag_incidence(3,indx_r);
            data(i).eigval = 1 ./ diag(jacob(indx_r, indx_c));
            data(i).rank = sum(abs(data(i).eigval) > 0);
            full_rank = (rcond(jacob(indx_r, indx_c)) > 1e-9);
        else
            data(i).eigval = [];
            data(i).rank = 0;
            full_rank = 1;
        end
        dr.eigval = [dr.eigval ; data(i).eigval];
        dr.full_rank = dr.full_rank && full_rank;
        %First order approximation
        if task ~= 1
            if (maximum_lag > 0)
                indx_r = find(M_.block_structure.block(i).lead_lag_incidence(3,:));
                indx_c = M_.block_structure.block(i).lead_lag_incidence(3,indx_r);
                ghx = - inv(jacob(indx_r, indx_c));
            end;
            ghu =  - inv(jacob(indx_r, indx_c)) * data(i).g1_x;
        end
      case 3
      %% ------------------------------------------------------------------
        %Solve Forward single equation
        if maximum_lag > 0 && n_pred > 0
            data(i).eigval = - jacob(1 , 1 : n_pred) / jacob(1 , n_pred + n_static + 1 : n_pred + n_static + n_pred + n_both);
            data(i).rank = 0;
        else
            data(i).eigval = [];
            data(i).rank = 0;
        end;
        dr.eigval = [dr.eigval ; data(i).eigval];
        %First order approximation
        if task ~= 1
            if (maximum_lag > 0)
                 ghx = - jacob(1 , 1 : n_pred) / jacob(1 , n_pred + n_static + 1 : n_pred + n_static + n_pred + n_both);
            else
                 ghx = 0;
            end;
            if other_endo_nbr
                fx = data(i).g1_o;
                % retrieves the derivatives with respect to endogenous
                % variable belonging to previous blocks
                fx_tm1 = zeros(n,other_endo_nbr);
                fx_t = zeros(n,other_endo_nbr);
                fx_tp1 = zeros(n,other_endo_nbr);
                % stores in fx_tm1 the lagged values of fx
                [r, c, lag] = find(data(i).lead_lag_incidence_other(1,:));
                fx_tm1(:,c) = fx(:,lag);
                % stores in fx the current values of fx
                [r, c, curr] = find(data(i).lead_lag_incidence_other(2,:));
                fx_t(:,c) = fx(:,curr);
                % stores in fx_tm1 the leaded values of fx
                [r, c, lead] = find(data(i).lead_lag_incidence_other(3,:));
                fx_tp1(:,c) = fx(:,lead);

                l_x = dr.ghx(data(i).other_endogenous,:);
                l_x_sv = dr.ghx(dr.state_var, 1:n_sv);
                
                selector_tm1 = M_.block_structure.block(i).tm1;
                ghx_other = - (fx_t * l_x + (fx_tp1 * l_x * l_x_sv) + fx_tm1 * selector_tm1) / jacob(1 , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                dr.ghx(endo, :) = dr.ghx(endo, :) + ghx_other;

            end;
            if exo_nbr
                fu = data(i).g1_x;
                if other_endo_nbr > 0
                    l_u_sv = dr.ghu(dr.state_var,:);
                    l_x = dr.ghx(data(i).other_endogenous,:);
                    l_u = dr.ghu(data(i).other_endogenous,:);
                    fu_complet = zeros(n, M_.exo_nbr);
                    fu_complet(:,data(i).exogenous) = fu;
                    ghu = -(fu_complet + fx_tp1 * l_x * l_u_sv + (fx_t) * l_u ) / jacob(1 , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                    exo = dr.exo_var;
                else
                    ghu = - fu  / jacob(1 , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                end;
            else
                 if other_endo_nbr > 0
                     l_u_sv = dr.ghu(dr.state_var,:);
                     l_x = dr.ghx(data(i).other_endogenous,:);
                     l_u = dr.ghu(data(i).other_endogenous,:);
                     ghu = -(fx_tp1 * l_x * l_u_sv + (fx_t) * l_u ) / jacob(1 , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                     exo = dr.exo_var;
                 else
                     ghu = [];
                 end
            end
        end
      case 4
      %% ------------------------------------------------------------------
        %Solve Backward single equation
        if maximum_lead > 0 && n_fwrd > 0
            data(i).eigval = - jacob(1 , n_pred + n - n_fwrd + 1 : n_pred + n) / jacob(1 , n_pred + n + 1 : n_pred + n + n_fwrd) ;
            data(i).rank = sum(abs(data(i).eigval) > 0);
            full_rank = (abs(jacob(1,n_pred+n+1: n_pred_n+n_fwrd)) > 1e-9);
        else
            data(i).eigval = [];
            data(i).rank = 0;
            full_rank = 1;
        end;
        dr.full_rank = dr.full_rank && full_rank;
        dr.eigval = [dr.eigval ; data(i).eigval];
      case 6
      %% ------------------------------------------------------------------
        %Solve Forward complete
        if maximum_lag > 0 && n_pred > 0
            data(i).eigval = eig(- jacob(: , 1 : n_pred) / ...
                                 jacob(: , (n_pred + n_static + 1 : n_pred + n_static + n_pred )));
            data(i).rank = 0;
            full_rank = (rcond(jacob(: , (n_pred + n_static + 1 : n_pred ...
                                          + n_static + n_pred ))) > 1e-9);
        else
            data(i).eigval = [];
            data(i).rank = 0;
            full_rank = 1;
        end;
        dr.eigval = [dr.eigval ; data(i).eigval];
        dr.full_rank = dr.full_rank && full_rank;
        if task ~= 1
            if (maximum_lag > 0)
                 ghx = - jacob(: , 1 : n_pred) / jacob(: , n_pred + n_static + 1 : n_pred + n_static + n_pred + n_both);
            else
                 ghx = 0;
            end;
            if other_endo_nbr
                fx = data(i).g1_o;
                % retrieves the derivatives with respect to endogenous
                % variable belonging to previous blocks
                fx_tm1 = zeros(n,other_endo_nbr);
                fx_t = zeros(n,other_endo_nbr);
                fx_tp1 = zeros(n,other_endo_nbr);
                % stores in fx_tm1 the lagged values of fx
                [r, c, lag] = find(data(i).lead_lag_incidence_other(1,:));
                fx_tm1(:,c) = fx(:,lag);
                % stores in fx the current values of fx
                [r, c, curr] = find(data(i).lead_lag_incidence_other(2,:));
                fx_t(:,c) = fx(:,curr);
                % stores in fx_tm1 the leaded values of fx
                [r, c, lead] = find(data(i).lead_lag_incidence_other(3,:));
                fx_tp1(:,c) = fx(:,lead);

                l_x = dr.ghx(data(i).other_endogenous,:);
                l_x_sv = dr.ghx(dr.state_var, 1:n_sv);
                
                selector_tm1 = M_.block_structure.block(i).tm1;
                ghx_other = - (fx_t * l_x + (fx_tp1 * l_x * l_x_sv) + fx_tm1 * selector_tm1) / jacob(: , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                dr.ghx(endo, :) = dr.ghx(endo, :) + ghx_other;
            end;
            if exo_nbr
                fu = data(i).g1_x;
                if other_endo_nbr > 0
                    l_u_sv = dr.ghu(dr.state_var,:);
                    l_x = dr.ghx(data(i).other_endogenous,:);
                    l_u = dr.ghu(data(i).other_endogenous,:);
                    fu_complet = zeros(n, M_.exo_nbr);
                    fu_complet(:,data(i).exogenous) = fu;
                    ghu = -(fu_complet + fx_tp1 * l_x * l_u_sv + (fx_t) * l_u ) / jacob(: , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                    exo = dr.exo_var;
                else
                    ghu = - fu  / jacob(: , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                end;
            else
                 if other_endo_nbr > 0
                     l_u_sv = dr.ghu(dr.state_var,:);
                     l_x = dr.ghx(data(i).other_endogenous,:);
                     l_u = dr.ghu(data(i).other_endogenous,:);
                     ghu = -(fx_tp1 * l_x * l_u_sv + (fx_t) * l_u ) / jacob(1 , n_pred + 1 : n_pred + n_static + n_pred + n_both);
                     exo = dr.exo_var;
                 else
                     ghu = [];
                 end
            end
        end
      case 7
      %% ------------------------------------------------------------------
        %Solve Backward complete
        if maximum_lead > 0 && n_fwrd > 0
            data(i).eigval = eig(- jacob(: , n_pred + n - n_fwrd + 1: n_pred + n))/ ...
                jacob(: , n_pred + n + 1 : n_pred + n + n_fwrd);
            data(i).rank = sum(abs(data(i).eigval) > 0);
            full_rank = (rcond(jacob(: , n_pred + n + 1 : n_pred + n + ...
                                     n_fwrd)) > 1e-9);
        else
            data(i).eigval = [];
            data(i).rank = 0;
            full_rank = 1;
        end;
        dr.full_rank = dr.full_rank && full_rank;
        dr.eigval = [dr.eigval ; data(i).eigval];
      case {5,8}
      %% ------------------------------------------------------------------
        %The lead_lag_incidence contains columns in the following order:
        %  static variables, backward variable, mixed variables and forward variables
        %  
        %Proceeds to a QR decomposition on the jacobian matrix in order to reduce the problem size
        index_c = lead_lag_incidence(2,:);             % Index of all endogenous variables present at time=t
        index_s = lead_lag_incidence(2,1:n_static);    % Index of all static endogenous variables present at time=t
        if n_static > 0
            [Q, junk] = qr(jacob(:,index_s));
            aa = Q'*jacob;
        else
            aa = jacob;
        end;
        index_0m = (n_static+1:n_static+n_pred) + indexi_0 - 1;
        index_0p = (n_static+n_pred+1:n) + indexi_0 - 1;
        index_m = 1:(n_pred+n_both);
        index_p  = lead_lag_incidence(3,find(lead_lag_incidence(3,:)));
        nyf = n_fwrd + n_both;
        A = aa(:,index_m);  % Jacobain matrix for lagged endogeneous variables
        B = aa(:,index_c);  % Jacobian matrix for contemporaneous endogeneous variables
        C = aa(:,index_p);  % Jacobain matrix for led endogeneous variables

        row_indx = n_static+1:n;

        if task ~= 1 && options_.dr_cycle_reduction == 1
            A1 = [aa(row_indx,index_m ) zeros(n_dynamic,n_fwrd)];
            B1 = [aa(row_indx,index_0m) aa(row_indx,index_0p) ];
            C1 = [zeros(n_dynamic,n_pred) aa(row_indx,index_p)];
            [ghx, info] = cycle_reduction(A1, B1, C1, options_.dr_cycle_reduction_tol);
            %ghx
            ghx = ghx(:,index_m);
            hx = ghx(1:n_pred+n_both,:);
            gx = ghx(1+n_pred:end,:);
        end
        
        if (task ~= 1 && ((options_.dr_cycle_reduction == 1 && info ==1) || options_.dr_cycle_reduction == 0)) || task == 1
            D = [[aa(row_indx,index_0m) zeros(n_dynamic,n_both) aa(row_indx,index_p)] ; [zeros(n_both, n_pred) eye(n_both) zeros(n_both, n_both + n_fwrd)]];
            E = [-aa(row_indx,[index_m index_0p])  ; [zeros(n_both, n_both + n_pred) eye(n_both, n_both + n_fwrd) ] ];
            D
            E

            [err, ss, tt, w, sdim, data(i).eigval, info1] = mjdgges(E,D,options_.qz_criterium);

            if (verbose)
                disp('eigval');
                disp(data(i).eigval);
            end;
            if info1
                info(1) = 2;
                info(2) = info1;
                return
            end
            nba = nd-sdim;
            if task == 1
                data(i).rank = rank(w(nd-nyf+1:end,nd-nyf+1:end));
                dr.full_rank = dr.full_rank && (rcond(w(nd-nyf+1:end,nd- ...
                                                        nyf+1:end)) > 1e-9);
                dr.eigval = [dr.eigval ; data(i).eigval];
            end
            if (verbose)
                disp(['sum eigval > 1 = ' int2str(sum(abs(data(i).eigval) > 1.)) ' nyf=' int2str(nyf) ' and dr.rank=' int2str(data(i).rank)]);
                disp(['data(' int2str(i) ').eigval']);
                disp(data(i).eigval);
            end;

            %First order approximation
            if task ~= 1
                if nba ~= nyf
                    sorted_roots = sort(abs(dr.eigval));
                    if isfield(options_,'indeterminacy_continuity')
                        if options_.indeterminacy_msv == 1
                            [ss,tt,w,q] = qz(e',d');
                            [ss,tt,w,junk] = reorder(ss,tt,w,q);
                            ss = ss';
                            tt = tt';
                            w  = w';
                            %nba = nyf;
                        end
                    else
                        if nba > nyf
                            temp = sorted_roots(nd-nba+1:nd-nyf)-1-options_.qz_criterium;
                            info(1) = 3;
                        elseif nba < nyf;
                            temp = sorted_roots(nd-nyf+1:nd-nba)-1-options_.qz_criterium;
                            info(1) = 4;
                        end
                        info(2) = temp'*temp;
                        return
                    end
                end
                indx_stable_root = 1: (nd - nyf);     %=> index of stable roots
                indx_explosive_root = n_pred + n_both + 1:nd;  %=> index of explosive roots
                % derivatives with respect to dynamic state variables
                % forward variables
                Z = w';
                Z11t = Z(indx_stable_root,    indx_stable_root)';
                Z21  = Z(indx_explosive_root, indx_stable_root);
                Z22  = Z(indx_explosive_root, indx_explosive_root);
                if ~isfloat(Z21) && (condest(Z21) > 1e9)
                    % condest() fails on a scalar under Octave
                    info(1) = 5;
                    info(2) = condest(Z21);
                    return;
                else
                    %gx = -inv(Z22) * Z21;
                    gx = - Z22 \ Z21;
                end

                % predetermined variables
                hx =  Z11t * inv(tt(indx_stable_root, indx_stable_root)) * ss(indx_stable_root, indx_stable_root) * inv(Z11t);

                k1 = 1:(n_pred+n_both);
                k2 = 1:(n_fwrd+n_both);

                ghx = [hx(k1,:); gx(k2(n_both+1:end),:)];
            end;
        end;
        
        if  task~= 1 
            %lead variables actually present in the model
            j4 = n_static+n_pred+1:n_static+n_pred+n_both+n_fwrd;   % Index on the forward and both variables
            j3 = nonzeros(lead_lag_incidence(2,j4)) - n_static - 2 * n_pred - n_both;  % Index on the non-zeros forward and both variables
            j4 = find(lead_lag_incidence(2,j4)); 
            
            if n_static > 0
                B_static = B(:,1:n_static);  % submatrix containing the derivatives w.r. to static variables
            else
                B_static = [];
            end;
            %static variables, backward variable, mixed variables and forward variables
            B_pred = B(:,n_static+1:n_static+n_pred+n_both);
            B_fyd = B(:,n_static+n_pred+n_both+1:end);

            % static variables
            if n_static > 0
                temp = - C(1:n_static,j3)*gx(j4,:)*hx;
                j5 = index_m;
                b = aa(:,index_c);
                b10 = b(1:n_static, 1:n_static);
                b11 = b(1:n_static, n_static+1:n);
                temp(:,j5) = temp(:,j5)-A(1:n_static,:);
                temp = b10\(temp-b11*ghx);
                ghx = [temp; ghx];
                temp = [];
            end;
            
            A_ = real([B_static C(:,j3)*gx+B_pred B_fyd]); % The state_variable of the block are located at [B_pred B_both]
            
            if other_endo_nbr
                if n_static > 0
                     fx = Q' * data(i).g1_o;
                else
                    fx = data(i).g1_o;
                end;
                % retrieves the derivatives with respect to endogenous
                % variable belonging to previous blocks
                fx_tm1 = zeros(n,other_endo_nbr);
                fx_t = zeros(n,other_endo_nbr);
                fx_tp1 = zeros(n,other_endo_nbr);
                % stores in fx_tm1 the lagged values of fx
                [r, c, lag] = find(data(i).lead_lag_incidence_other(1,:));
                fx_tm1(:,c) = fx(:,lag);
                % stores in fx the current values of fx
                [r, c, curr] = find(data(i).lead_lag_incidence_other(2,:));
                fx_t(:,c) = fx(:,curr);
                % stores in fx_tp1 the leaded values of fx
                [r, c, lead] = find(data(i).lead_lag_incidence_other(3,:));
                fx_tp1(:,c) = fx(:,lead);

                l_x = dr.ghx(data(i).other_endogenous,:);
                
                l_x_sv = dr.ghx(dr.state_var, :);
                
                selector_tm1 = M_.block_structure.block(i).tm1; 

                B_ = [zeros(size(B_static)) zeros(n,n_pred) C(:,j3) ];
                C_ = l_x_sv;
                D_ = (fx_t * l_x + fx_tp1 * l_x * l_x_sv + fx_tm1 * selector_tm1 );
                % Solve the Sylvester equation:
                % A_ * gx + B_ * gx * C_ + D_ = 0
                if block_type == 5
                    vghx_other = - inv(kron(eye(size(D_,2)), A_) + kron(C_', B_)) * vec(D_);
                    ghx_other = reshape(vghx_other, size(D_,1), size(D_,2));
                elseif options_.sylvester_fp == 1
                    ghx_other = gensylv_fp(A_, B_, C_, D_, i, options_.sylvester_fixed_point_tol);
                else 
                    [err, ghx_other] = gensylv(1, A_, B_, C_, -D_);
                end;
                if options_.aim_solver ~= 1 && options_.use_qzdiv
                   % Necessary when using Sims' routines for QZ
                   ghx_other = real(ghx_other);
                end
                
                dr.ghx(endo, :) = dr.ghx(endo, :) + ghx_other;
            end;

            if exo_nbr
                if n_static > 0
                    fu = Q' * data(i).g1_x;
                else
                    fu = data(i).g1_x;
                end;

                if other_endo_nbr > 0
                    l_u_sv = dr.ghu(dr.state_var,:);
                    l_x = dr.ghx(data(i).other_endogenous,:);
                    l_u = dr.ghu(data(i).other_endogenous,:);
                    fu_complet = zeros(n, M_.exo_nbr);
                    fu_complet(:,data(i).exogenous) = fu;
                    % Solve the equation in ghu:
                    % A_ * ghu + (fu_complet + fx_tp1 * l_x * l_u_sv + (fx_t + B_ * ghx_other) * l_u ) = 0
                    
                    ghu = -A_\ (fu_complet + fx_tp1 * l_x * l_u_sv + fx_t * l_u + B_ * ghx_other  * l_u_sv  );
                    exo = dr.exo_var;
                else
                    ghu = - A_ \ fu;
                end;
            else
                if other_endo_nbr > 0
                    l_u_sv = dr.ghu(dr.state_var,:);
                    l_x = dr.ghx(data(i).other_endogenous,:);
                    l_u = dr.ghu(data(i).other_endogenous,:);
                    % Solve the equation in ghu:
                    % A_ * ghu + (fx_tp1 * l_x * l_u_sv + (fx_t + B_ * ghx_other) * l_u ) = 0
                    ghu = -real(A_)\ (fx_tp1 * l_x * l_u_sv + (fx_t * l_u + B_ * ghx_other * l_u_sv) );
                    exo = dr.exo_var;
                else
                    ghu = [];
                end;
            end


            
            if options_.loglinear == 1
                error('log linear option is for the moment not supported in first order approximation for a block decomposed mode');
%                 k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
%                 klag = dr.kstate(k,[1 2]);
%                 k1 = dr.order_var;
%                 
%                 ghx = repmat(1./dr.ys(k1),1,size(ghx,2)).*ghx.* ...
%                       repmat(dr.ys(k1(klag(:,1)))',size(ghx,1),1);
%                 ghu = repmat(1./dr.ys(k1),1,size(ghu,2)).*ghu;
            end


            if options_.aim_solver ~= 1 && options_.use_qzdiv
                % Necessary when using Sims' routines for QZ
                ghx = real(ghx);
                ghu = real(ghu);
            end

            %exogenous deterministic variables
            if exo_det_nbr > 0
                error('deterministic exogenous are not yet implemented in first order approximation for a block decomposed model');
%                 f1 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+2:end,order_var))));
%                 f0 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var))));
%                 fudet = data(i).g1_xd;
%                 M1 = inv(f0+[zeros(n,n_static) f1*gx zeros(n,nyf-n_both)]);
%                 M2 = M1*f1;
%                 dr.ghud = cell(M_.exo_det_length,1);
%                 dr.ghud{1} = -M1*fudet;
%                 for i = 2:M_.exo_det_length
%                     dr.ghud{i} = -M2*dr.ghud{i-1}(end-nyf+1:end,:);
%                 end
            end
        end
    end;
    if task ~=1
        if (maximum_lag > 0 && (n_pred > 0 || n_both > 0))
            sorted_col_dr_ghx = M_.block_structure.block(i).sorted_col_dr_ghx;
            dr.ghx(endo, sorted_col_dr_ghx) = dr.ghx(endo, sorted_col_dr_ghx) + ghx;
            data(i).ghx = ghx;
            data(i).pol.i_ghx = sorted_col_dr_ghx;
        else
            data(i).pol.i_ghx = [];
        end;
        data(i).ghu = ghu;
        dr.ghu(endo, exo) = ghu;
        data(i).pol.i_ghu = exo;
    end;
    
   if (verbose)
        disp('dr.ghx');
        dr.ghx
        disp('dr.ghu');
        dr.ghu
   end; 
   
end;
M_.block_structure.block = data ;
if (verbose)
        disp('dr.ghx');
        disp(real(dr.ghx));
        disp('dr.ghu');
        disp(real(dr.ghu));
end; 
if (task == 1)
    return;
end;
