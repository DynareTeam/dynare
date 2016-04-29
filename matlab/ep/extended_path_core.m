function [y, info_convergence, endogenousvariablespaths] = extended_path_core(periods,endo_nbr,exo_nbr,positive_var_indx, ...
                                exo_simul,init,initial_conditions,...
                                steady_state, ...
                                debug,bytecode_flag,order,M,pfm,algo,solve_algo,stack_solve_algo,...
                                olmmcp,options,oo,initialguess)

% Copyright (C) 2016 Dynare Team
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

ep = options.ep;
if init% Compute first order solution (Perturbation)...
    endo_simul = simult_(initial_conditions,oo.dr,exo_simul(2:end,:),1);
else
    if nargin==20 && ~isempty(initialguess)
        % Note that the first column of initialguess should be equal to initial_conditions.
        endo_simul = initialguess;
    else
        endo_simul = [initial_conditions repmat(steady_state,1,periods+1)];
    end
end

oo.endo_simul = endo_simul;

if debug
    save ep_test_1 endo_simul exo_simul
end

if bytecode_flag && ~ep.stochastic.order
    [flag,tmp] = bytecode('dynamic',endo_simul,exo_simul, M_.params, endo_simul, periods);
else
    flag = 1;
end

if flag
    if order == 0
        options.periods = periods;
        options.block = pfm.block;
        oo.endo_simul = endo_simul;
        oo.exo_simul = exo_simul;
        oo.steady_state = steady_state;
        options.bytecode = bytecode_flag;
        options.lmmcp = olmmcp;
        options.solve_algo = solve_algo;
        options.stack_solve_algo = stack_solve_algo;
        [tmp, maxerror] = perfect_foresight_solver_core(M, options, oo);
        if maxerror>options.dynatol.f
            flag = false;
        else
            flag = true;
        end
        if ~flag && ~options.no_homotopy
            exo_orig = oo.exo_simul;
            endo_simul = repmat(steady_state,1,periods+1);
            for i = 1:10
                weight = i/10;
                oo.endo_simul = [weight*initial_conditions + (1-weight)*steady_state ...
                                 endo_simul];
                oo.exo_simul = repmat((1-weight)*oo.exo_steady_state', ...
                                      size(oo.exo_simul,1),1) + weight*exo_orig;
                [tmp,flag] = perfect_foresight_solver_core(M,options,oo);
                disp([i,flag])
                if ~flag
                    break
                end
                endo_simul = tmp.endo_simul;
            end
        end
        info_convergence = flag;
    else
        switch(algo)
          case 0
            [flag,endo_simul] = ...
                solve_stochastic_perfect_foresight_model(endo_simul,exo_simul,pfm,ep.stochastic.quadrature.nodes,ep.stochastic.order);
          case 1
            [flag,endo_simul] = ...
                solve_stochastic_perfect_foresight_model_1(endo_simul,exo_simul,options,pfm,ep.stochastic.order);
        end
        tmp.endo_simul = endo_simul;
        info_convergence = ~flag;
    end
end
if info_convergence
    y = tmp.endo_simul(:,2);
else
    y = NaN(size(endo_nbr,1));
end
