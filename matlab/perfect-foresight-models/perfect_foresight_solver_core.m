function [oo_, maxerror] = perfect_foresight_solver_core(M_, options_, oo_)
%function [oo_, maxerror] = perfect_foresight_solver_core(M_, options_, oo_)
% Core function calling solvers for perfect foresight model
%
% INPUTS 
% - M_                  [struct] contains a description of the model.
% - options_            [struct] contains various options.
% - oo_                 [struct] contains results
%
% OUTPUTS 
% - oo_                 [struct] contains results
% - maxerror            [double] contains the maximum absolute error

% Copyright (C) 2015-2016 Dynare Team
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

if options_.lmmcp.status
    options_.stack_solve_algo=7;
    options_.solve_algo = 10;
end

if options_.linear_approximation && ~(isequal(options_.stack_solve_algo,0) || isequal(options_.stack_solve_algo,7))
    error('perfect_foresight_solver: Option linear_approximation is only available with option stack_solve_algo equal to 0.')
end

if options_.linear && isequal(options_.stack_solve_algo,0)
    options_.linear_approximation = 1;
end

if options_.block
    if options_.bytecode
        try
            [info, tmp] = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state,1,options_.periods+2), options_.periods);
        catch
            info = 1;
        end
        if info
            oo_.deterministic_simulation.status = false;
        else
            oo_.endo_simul = tmp;
            oo_.deterministic_simulation.status = true;
        end
        if options_.no_homotopy
            mexErrCheck('bytecode', info);
        end
    else
        oo_ = feval([M_.fname '_dynamic'], options_, M_, oo_);
    end
else
    if options_.bytecode
        try
            [info, tmp] = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state,1,options_.periods+2), options_.periods);
        catch
            info = 1;
        end
        if info
            oo_.deterministic_simulation.status = false;
        else
            oo_.endo_simul = tmp;
            oo_.deterministic_simulation.status = true;
        end
        if options_.no_homotopy
            mexErrCheck('bytecode', info);
        end
    else
        if M_.maximum_endo_lead == 0 % Purely backward model
            [oo_.endo_simul, oo_.deterministic_simulation] = ...
                sim1_purely_backward(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 % Purely forward model
            [oo_.endo_simul, oo_.deterministic_simulation] = ...
                sim1_purely_forward(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        else % General case
            switch options_.stack_solve_algo
              case 0
                if options_.linear_approximation
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        sim1_linear(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
                else
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
                end
              case 6
                if options_.linear_approximation
                    error('Invalid value of stack_solve_algo option!')
                end
                [oo_.endo_simul, oo_.deterministic_simulation] = ...
                    sim1_lbj(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
              case 7
                if options_.linear_approximation
                    if isequal(options_.solve_algo, 10)
                        warning('It would be more efficient to set option solve_algo equal to 0!')
                    end
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        solve_stacked_linear_problem(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
                else
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        solve_stacked_problem(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
                end
              otherwise
                error('Invalid value of stack_solve_algo option!')
            end
        end
    end
end

if nargout>1
    y0 = oo_.endo_simul(:,1);
    yT = oo_.endo_simul(:,options_.periods+2);
    yy  = oo_.endo_simul(:,2:options_.periods+1);
    if ~exist('illi')
        illi = M_.lead_lag_incidence';
        [i_cols,junk,i_cols_j] = find(illi(:));
        illi = illi(:,2:3);
        [i_cols_J1,junk,i_cols_1] = find(illi(:));
        i_cols_T = nonzeros(M_.lead_lag_incidence(1:2,:)');
    end
    if options_.block && ~options_.bytecode
        maxerror = oo_.deterministic_simulation.error;
    else
        if options_.bytecode
            [chck, residuals, junk]= bytecode('dynamic','evaluate', oo_.endo_simul, oo_.exo_simul, M_.params, oo_.steady_state, 1);
        else
            residuals = perfect_foresight_problem(yy(:),str2func([M_.fname '_dynamic']), y0, yT, ...
                                                  oo_.exo_simul,M_.params,oo_.steady_state, ...
                                                  M_.maximum_lag,options_.periods,M_.endo_nbr,i_cols, ...
                                                  i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                                                  M_.NNZDerivatives(1));
        end
        maxerror = max(max(abs(residuals)));
    end
end


