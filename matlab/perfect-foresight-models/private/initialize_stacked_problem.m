function [options, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, ...
          dynamicmodel] = initialize_stacked_problem(endogenousvariables, options, M, steadystate_y)
% function [options, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, ...
%           dynamicmodel] = initialize_stacked_problem(endogenousvariables, options, M, steadystate_y)
% Sets up the stacked perfect foresight problem for use with dynare_solve.m
%
% INPUTS 
% - endogenousvariables [double] N*T array, paths for the endogenous variables (initial guess).
% - options             [struct] contains various options.
% - M                   [struct] contains a description of the model.
% - steadystate_y       [double] N*1 array, steady state for the endogenous variables.
% OUTPUTS 
% - options             [struct] contains various options.
% - y0                  [double] N*1 array, initial conditions for the endogenous variables
% - yT                  [double] N*1 array, terminal conditions for the endogenous variables
% - z                   [double] T*M array, paths for the exogenous variables.
% - i_cols              [double] indices of variables appearing in M.lead_lag_incidence
%                                and that need to be passed to _dynamic-file
% - i_cols_J1           [double] indices of contemporaneous and forward looking variables 
%                                appearing in M.lead_lag_incidence
% - i_cols_T            [double] columns of dynamic Jacobian related to
%                                contemporaneous and backward-looking
%                                variables (relevant in last period)
% - i_cols_j            [double] indices of variables in M.lead_lag_incidence
%                                in dynamic Jacobian (relevant in intermediate periods)
% - i_cols_1            [double] indices of contemporaneous and forward looking variables in
%                                M.lead_lag_incidence in dynamic Jacobian (relevant in first period)
% - dynamicmodel        [handle] function handle to _dynamic-file

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
    
periods = options.periods;
if (options.solve_algo == 10)
    if ~isfield(options.lmmcp,'lb')
        [lb,ub,pfm.eq_index] = get_complementarity_conditions(M,options.ramsey_policy);
        if options.linear_approximation
            lb = lb - steadystate_y;
            ub = ub - steadystate_y;
        end
        options.lmmcp.lb = repmat(lb,periods,1);
        options.lmmcp.ub = repmat(ub,periods,1);
    end
elseif (options.solve_algo == 11)
    if ~isfield(options.mcppath,'lb')
        [lb,ub,pfm.eq_index] = get_complementarity_conditions(M,options.ramsey_policy);
        if options.linear_approximation
            lb = lb - steadystate_y;
            ub = ub - steadystate_y;
        end
        options.mcppath.lb = repmat(lb,periods,1);
        options.mcppath.ub = repmat(ub,periods,1);
    end
end

y0 = endogenousvariables(:,1);
yT = endogenousvariables(:,periods+2);
z = endogenousvariables(:,2:periods+1);
illi = M.lead_lag_incidence';
[i_cols, junk,i_cols_j] = find(illi(:));
illi = illi(:,2:3);
[i_cols_J1, junk,i_cols_1] = find(illi(:));
i_cols_T = nonzeros(M.lead_lag_incidence(1:2,:)');
dynamicmodel = str2func([M.fname,'_dynamic']);