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
%   ys                        vector           steady state
%   params                    vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2013 Dynare Team
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
        ys_init = h_set_auxiliary_variables(ys_init,exo_ss,M.params);
    end

    if options.ramsey_policy
        [ys,params] = dyn_ramsey_static(ys_init,M,options,oo);
    elseif steadystate_flag
        % explicit steady state file
        [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M, ...
                                                       options);
        if info(1)
            return;
        end
    elseif (options.bytecode == 0 && options.block == 0)
        if options.linear == 0
            % non linear model
            [ys,check] = dynare_solve([M.fname '_static'],...
                                      ys_init,...
                                      options.jacobian_flag, ...
                                      exo_ss, params);
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
                disp('Check whether your model in truly linear')
            elseif isempty(ii) && max(abs(fvec)) > 1e-12
                ys = ys_init-jacob\fvec;
            else
                ys = ys_init;
            end

        end
    else
        % block or bytecode
        [ys,check] = dynare_solve_block_or_bytecode(ys_init,exo_ss, params, ...
                                                    options, M);
    end

    if check
        if options.steadystate_flag
            info(1)= 19;
            resid = check1 ;
        else
            info(1)= 20;
            resid = evaluate_static_model(ys_init,exo_ss,params,M,options);
        end
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
        zx = repmat([oo.exo_simul oo.exo_det_simul],M.maximum_lead + M.maximum_lag + 1, 1);
        if options.bytecode
            [chck, r, junk]= bytecode('dynamic','evaluate', z, zx, M.params, ys, 1);
            mexErrCheck('bytecode', chck);
        elseif options.block
            [r, data] = feval([M.fname '_dynamic'], z', zx, M.params, ys, M.maximum_lag+1, data);
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

