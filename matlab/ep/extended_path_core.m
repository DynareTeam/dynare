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
    save ep_test_1.mat endo_simul exo_simul
end

if bytecode_flag && ~ep.stochastic.order
    [flag, tmp] = bytecode('dynamic', endo_simul, exo_simul, M_.params, endo_simul, periods);
else
    flag = true;
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
        tmp = perfect_foresight_solver_core(M, options, oo);
        if ~tmp.deterministic_simulation.status 
            info_convergence = false;
        else
            info_convergence = true;
        end
    else
        switch(algo)
          case 0
            [flag, tmp.endo_simul] = ...
                solve_stochastic_perfect_foresight_model(endo_simul, exo_simul, pfm, ep.stochastic.quadrature.nodes, ep.stochastic.order);
          case 1
            [flag, tmp.endo_simul] = ...
                solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, options, pfm, ep.stochastic.order);
        end
        info_convergence = ~flag;
    end
end

if ~info_convergence && ~options.no_homotopy
    [info_convergence, tmp.endo_simul] = extended_path_homotopy(endo_simul, exo_simul, M, options, oo, pfm, ep, order, algo, 2, debug);
end

if info_convergence
    y = tmp.endo_simul(:,2);
else
    y = NaN(size(endo_nbr,1));
end

endogenousvariablespaths = tmp.endo_simul;