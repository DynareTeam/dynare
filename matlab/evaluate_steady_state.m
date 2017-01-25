function [ys,params,info] = evaluate_steady_state(ys_init,M,options,oo,steadystate_check_flag)
% function [ys,params,info] = evaluate_steady_state(ys_init,M,options,oo,steadystate_check_flag)
% Computes the steady state
%
% INPUTS
%   ys_init                   vector           initial values used to compute the steady
%                                                 state
%   M                         struct           model structure
%   options                   struct           options
%   oo                        struct           output results
%   steadystate_check_flag    boolean          if true, check that the
%                                              steadystate verifies the
%                                              static model
%
% OUTPUTS
%   ys                        vector           steady state (in declaration order)
%   params                    vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2016 Dynare Team
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
    check = 0;

    steadystate_flag = options.steadystate_flag;
    params = M.params;
    exo_ss = [oo.exo_steady_state; oo.exo_det_steady_state];

    if length(M.aux_vars) > 0
        h_set_auxiliary_variables = str2func([M.fname '_set_auxiliary_variables']);
        if ~steadystate_flag
            ys_init = h_set_auxiliary_variables(ys_init,exo_ss,M.params);
        end
    end

    if options.ramsey_policy
        if steadystate_flag
            % explicit steady state file
            [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M, ...
                                                           options,steadystate_check_flag);
            %test whether it solves model conditional on the instruments
            resids = evaluate_static_model(ys,exo_ss,params,M,options);
            n_multipliers=M.ramsey_eq_nbr;
            nan_indices=find(isnan(resids(n_multipliers+1:end)));

            if ~isempty(nan_indices)
                if options.debug
                    fprintf('\nevaluate_steady_state: The steady state file computation for the Ramsey problem resulted in NaNs.\n')
                    fprintf('evaluate_steady_state: The steady state was computed conditional on the following initial instrument values: \n')
                    for ii = 1:size(options.instruments,1);
                        fprintf('\t %s \t %f \n',options.instruments(ii,:),ys_init(strmatch(options.instruments(ii,:),M.endo_names,'exact')))
                    end
                    fprintf('evaluate_steady_state: The problem occured in the following equations: \n')
                    fprintf('\t Equation(s): ')
                    for ii=1:length(nan_indices)
                        fprintf('%d, ',nan_indices(ii));
                    end
                    skipline();
                    fprintf('evaluate_steady_state: If those initial values are not admissable, change them using an initval-block.\n')
                    skipline(2);
                end
                info(1) = 84;
                info(2) = resids'*resids;
                return;
            end
            
            if any(imag(ys(n_multipliers+1:end)))
                if options.debug
                    fprintf('\nevaluate_steady_state: The steady state file computation for the Ramsey problem resulted in complex numbers.\n')
                    fprintf('evaluate_steady_state: The steady state was computed conditional on the following initial instrument values: \n')
                    for ii = 1:size(options.instruments,1);
                        fprintf('\t %s \t %f \n',options.instruments(ii,:),ys_init(strmatch(options.instruments(ii,:),M.endo_names,'exact')))
                    end
                    fprintf('evaluate_steady_state: If those initial values are not admissable, change them using an initval-block.\n')
                    skipline(2);
                end
                info(1) = 86;
                info(2) = resids'*resids;
                return;
            end

            if max(abs(resids(n_multipliers+1:end))) > options.solve_tolf %does it solve for all variables except for the Lagrange multipliers
                if options.debug
                    fprintf('\nevaluate_steady_state: The steady state file does not solve the steady state for the Ramsey problem.\n')
                    fprintf('evaluate_steady_state: Conditional on the following instrument values: \n')
                    for ii = 1:size(options.instruments,1);
                        fprintf('\t %s \t %f \n',options.instruments(ii,:),ys_init(strmatch(options.instruments(ii,:),M.endo_names,'exact')))
                    end
                    fprintf('evaluate_steady_state: the following equations have non-zero residuals: \n')
                    for ii=n_multipliers+1:M.endo_nbr
                        if abs(resids(ii)) > options.solve_tolf
                            fprintf('\t Equation number %d: %f\n',ii-n_multipliers, resids(ii))
                        end
                    end
                    skipline(2);
                end
                info(1) = 85;
                info(2) = resids'*resids;
                return;
            end
        end
        if options.debug
            infrow=find(isinf(ys_init));
            if ~isempty(infrow)
                fprintf('\nevaluate_steady_state: The initial values for the steady state of the following variables are Inf:\n');
                for iter=1:length(infrow)
                    fprintf('%s\n',M.endo_names(infrow(iter),:));
                end
            end
            nanrow=find(isnan(ys_init));
            if ~isempty(nanrow)
                fprintf('\nevaluate_steady_state: The initial values for the steady state of the following variables are NaN:\n');
                for iter=1:length(nanrow)
                    fprintf('%s\n',M.endo_names(nanrow(iter),:));
                end
            end
        end
        %either if no steady state file or steady state file without problems
        [ys,params,info] = dyn_ramsey_static(ys_init,M,options,oo);
        if info
           return;
        end
        %check whether steady state really solves the model
        resids = evaluate_static_model(ys,exo_ss,params,M,options);

        n_multipliers=M.ramsey_eq_nbr;
        nan_indices_multiplier=find(isnan(resids(1:n_multipliers)));
        nan_indices=find(isnan(resids(n_multipliers+1:end)));

        if ~isempty(nan_indices)
            if options.debug
                fprintf('\nevaluate_steady_state: The steady state computation for the Ramsey problem resulted in NaNs.\n')
                fprintf('evaluate_steady_state: The steady state computation resulted in the following instrument values: \n')
                for i = 1:size(options.instruments,1);
                    fprintf('\t %s \t %f \n',options.instruments(i,:),ys(strmatch(options.instruments(i,:),M.endo_names,'exact')))
                end
                fprintf('evaluate_steady_state: The problem occured in the following equations: \n')
                fprintf('\t Equation(s): ')
                for ii=1:length(nan_indices)
                    fprintf('%d, ',nan_indices(ii));
                end
                skipline();
            end
            info(1) = 82;
            return;
        end

        if ~isempty(nan_indices_multiplier)
            if options.debug
                fprintf('\nevaluate_steady_state: The steady state computation for the Ramsey problem resulted in NaNs in the auxiliary equations.\n')
                fprintf('evaluate_steady_state: The steady state computation resulted in the following instrument values: \n')
                for i = 1:size(options.instruments,1);
                    fprintf('\t %s \t %f \n',options.instruments(i,:),ys(strmatch(options.instruments(i,:),M.endo_names,'exact')))
                end
                fprintf('evaluate_steady_state: The problem occured in the following equations: \n')
                fprintf('\t Auxiliary equation(s): ')
                for ii=1:length(nan_indices_multiplier)
                    fprintf('%d, ',nan_indices_multiplier(ii));
                end
                skipline();
            end
            info(1) = 83;
            return;
        end

        if max(abs(resids)) > options.solve_tolf %does it solve for all variables including the auxiliary ones
            if options.debug
                fprintf('\nevaluate_steady_state: The steady state for the Ramsey problem could not be computed.\n')
                fprintf('evaluate_steady_state: The steady state computation stopped with the following instrument values:: \n')
                for i = 1:size(options.instruments,1);
                    fprintf('\t %s \t %f \n',options.instruments(i,:),ys(strmatch(options.instruments(i,:),M.endo_names,'exact')))
                end
                fprintf('evaluate_steady_state: The following equations have non-zero residuals: \n')
                for ii=1:n_multipliers
                    if abs(resids(ii)) > options.solve_tolf/100
                        fprintf('\t Auxiliary Ramsey equation number %d: %f\n',ii, resids(ii))
                    end
                end
                for ii=n_multipliers+1:M.endo_nbr
                    if abs(resids(ii)) > options.solve_tolf/100
                        fprintf('\t Equation number %d: %f\n',ii-n_multipliers, resids(ii))
                    end
                end
                skipline(2);
            end
            info(1) = 81;
            info(2) = resids'*resids;
            return;
        end
    elseif steadystate_flag
        % explicit steady state file
        [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M, ...
                                                       options,steadystate_check_flag);
        if size(ys,2)>size(ys,1)
            error('STEADY: steady_state-file must return a column vector, not a row vector.')
        end
        if info(1)
            return;
        end
    elseif (options.bytecode == 0 && options.block == 0)
        if options.linear == 0
            % non linear model
            static_model = str2func([M.fname '_static']);
            [ys,check] = dynare_solve(@static_problem,...
                                      ys_init,...
                                      options, exo_ss, params,...
                                      M.endo_nbr,...
                                      static_model);
            if check && options.debug
                [ys,check,fvec,fjac] = dynare_solve(@static_problem,...
                          ys_init,...
                          options, exo_ss, params,...
                          M.endo_nbr,...
                          static_model);
                [infrow,infcol]=find(isinf(fjac) | isnan(fjac));
                if ~isempty(infrow)
                    fprintf('\nSTEADY:  The Jacobian at the initial values contains Inf or NaN. The problem arises from: \n')
                    display_problematic_vars_Jacobian(infrow,infcol,M,ys_init,'static','STEADY: ')
                end
                problematic_equation = find(~isfinite(fvec));                
                if ~isempty(problematic_equation)
                    fprintf('\nSTEADY:  numerical initial values or parameters incompatible with the following equations\n')
                    disp(problematic_equation')
                    fprintf('Please check for example\n')
                    fprintf('   i) if all parameters occurring in these equations are defined\n')
                    fprintf('  ii) that no division by an endogenous variable initialized to 0 occurs\n')
                end
            end
        else
            % linear model
            fh_static = str2func([M.fname '_static']);
            [fvec,jacob] = fh_static(ys_init,exo_ss, ...
                                     params);

            ii = find(~isfinite(fvec));
            if ~isempty(ii)
                ys=fvec;
                check=1;
                disp(['STEADY:  numerical initial values or parameters incompatible with the following' ...
                      ' equations'])
                disp(ii')
                disp('Check whether your model is truly linear. Put "resid(1);" before "steady;" to see the problematic equations.\n')
            elseif isempty(ii) && max(abs(fvec)) > 1e-12
                ys = ys_init-jacob\fvec;
                resid = evaluate_static_model(ys,exo_ss,params,M,options);
                if max(abs(resid)) > 1e-6
                    check=1;
                    fprintf('STEADY: No steady state for your model could be found\n')
                    fprintf('STEADY: Check whether your model is truly linear. Put "resid(1);" before "steady;" to see the problematic equations.\n')
                end

            else
                ys = ys_init;
            end
            if options.debug
                if any(any(isinf(jacob) | isnan(jacob)))
                    [infrow,infcol]=find(isinf(jacob) | isnan(jacob));
                    fprintf('\nSTEADY:  The Jacobian contains Inf or NaN. The problem arises from: \n\n')
                    for ii=1:length(infrow)
                        fprintf('STEADY:  Derivative of Equation %d with respect to Variable %s  (initial value of %s: %g) \n',infrow(ii),deblank(M.endo_names(infcol(ii),:)),deblank(M.endo_names(infcol(ii),:)),ys_init(infcol(ii)))
                    end
                    fprintf('Check whether your model is truly linear. Put "resid(1);" before "steady;" to see the problematic equations.\n')
                end
            end
        end
    else
        % block or bytecode
        [ys,check] = dynare_solve_block_or_bytecode(ys_init,exo_ss, params, ...
                                                    options, M);
    end

    if check
        info(1)= 20;
        %make sure ys contains auxiliary variables in case of problem with dynare_solve
        if length(M.aux_vars) > 0 && ~steadystate_flag
            ys = h_set_auxiliary_variables(ys,exo_ss,M.params);
        end
        resid = evaluate_static_model(ys,exo_ss,params,M,options);
        info(2) = resid'*resid ;
        if isnan(info(2))
            info(1)=22;
        end
        return
    end

    % If some equations are tagged [static] or [dynamic], verify consistency
    if M.static_and_dynamic_models_differ
        % Evaluate residual of *dynamic* model using the steady state
        % computed on the *static* one
        z = repmat(ys,1,M.maximum_lead + M.maximum_lag + 1);
        zx = repmat([exo_ss'], M.maximum_lead + M.maximum_lag + 1, 1);
        if options.bytecode
            [chck, r, junk]= bytecode('dynamic','evaluate', z, zx, M.params, ys, 1);
            mexErrCheck('bytecode', chck);
        elseif options.block
            [r, oo.dr] = feval([M.fname '_dynamic'], z', zx, M.params, ys, M.maximum_lag+1, oo.dr);
        else
            iyv = M.lead_lag_incidence';
            iyr0 = find(iyv(:));
            xys = z(iyr0);
            r = feval([M.fname '_dynamic'], z(iyr0), zx, M.params, ys, M.maximum_lag + 1);
        end

        % Fail if residual greater than tolerance
        if max(abs(r)) > options.solve_tolf
            info(1) = 25;
            return
        end
    end

    if ~isreal(ys)
        info(1) = 21;
        info(2) = sum(imag(ys).^2);
        ys = real(ys);
        return
    end

    if ~isempty(find(isnan(ys)))
        info(1) = 22;
        info(2) = NaN;
        return
    end

function [resids,jac] = static_problem(y,x,params,nvar,fh_static_model)
    [r,j] = fh_static_model(y,x,params);
    resids = r(1:nvar);
    jac = j(1:nvar,1:nvar);
