function [endogenousvariables, info] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, M, options);
% [endogenousvariables, info] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, M, options);
% Solves the perfect foresight model using dynare_solve
%
% INPUTS 
% - endogenousvariables [double] N*T array, paths for the endogenous variables (initial guess).
% - exogenousvariables  [double] T*M array, paths for the exogenous variables.
% - steadystate         [double] N*1 array, steady state for the endogenous variables.
% - M                   [struct] contains a description of the model.
% - options             [struct] contains various options.
%
% OUTPUTS 
% - endogenousvariables [double] N*T array, paths for the endogenous variables (solution of the perfect foresight model).
% - info                [struct] contains informations about the results.

% Copyright (C) 2015-16 Dynare Team
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
    
[options, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, dynamicmodel] = ...
    initialize_stacked_problem(endogenousvariables, options, M, steadystate);

if (options.solve_algo == 10 || options.solve_algo == 11)% mixed complementarity problem
    [lb,ub,eq_index] = get_complementarity_conditions(M,options.ramsey_policy);
    if options.linear_approximation
        lb = lb - steadystate_y;
        ub = ub - steadystate_y;
    end
    if options.solve_algo == 10
        options.lmmcp.lb = repmat(lb,options.periods,1);
        options.lmmcp.ub = repmat(ub,options.periods,1);
    elseif options.solve_algo == 11
        options.mcppath.lb = repmat(lb,options.periods,1);
        options.mcppath.ub = repmat(ub,options.periods,1);
    end
    [y, check] = dynare_solve(@perfect_foresight_mcp_problem,z(:),options, ...
                              dynamicmodel, y0, yT, ...
                              exogenousvariables, M.params, steadystate, ...
                              M.maximum_lag, options.periods, M.endo_nbr, i_cols, ...
                              i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                              M.NNZDerivatives(1),eq_index);
else
    [y, check] = dynare_solve(@perfect_foresight_problem,z(:),options, ...
                              dynamicmodel, y0, yT, ...
                              exogenousvariables, M.params, steadystate, ...
                              M.maximum_lag, options.periods, M.endo_nbr, i_cols, ...
                              i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                              M.NNZDerivatives(1));
end

if all(imag(y)<.1*options.dynatol.x)
    if ~isreal(y)
        y = real(y);
    end
else
    check = 1;
end

endogenousvariables = [y0 reshape(y, M.endo_nbr, options.periods) yT];

if check
    info.status = false;
else
    
    info.status = true;
end