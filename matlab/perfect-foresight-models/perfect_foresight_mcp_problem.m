function [residuals,JJacobian] = perfect_foresight_mcp_problem(y, dynamic_function, Y0, YT, ...
                                           exo_simul, params, steady_state, ...
                                           maximum_lag, T, ny, i_cols, ...
                                           i_cols_J1, i_cols_1, i_cols_T, ...
                                           i_cols_j,nnzJ,eq_index)
% function [residuals,JJacobian] = perfect_foresight_mcp_problem(x, model_dynamic, Y0, YT,exo_simul,
% params, steady_state, maximum_lag, periods, ny, i_cols,i_cols_J1, i_cols_1,
% i_cols_T, i_cols_j, nnzA) 
% computes the residuals and th Jacobian matrix
% for a perfect foresight problem over T periods.
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 1996-2015 Dynare Team
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


    YY = [Y0; y; YT];
    
    residuals = zeros(T*ny,1);
    if nargout == 2
        iJacobian = cell(T,1);
    end

    i_rows = 1:ny;
    offset = 0;
    i_cols_J = i_cols;

    for it = 2:(T+1)
        if nargout == 1
             res = dynamic_function(YY(i_cols),exo_simul, params, ...
                                                         steady_state,it);
             residuals(i_rows) = res(eq_index);
        elseif nargout == 2
            [res,jacobian] = dynamic_function(YY(i_cols),exo_simul, params, ...
                                                         steady_state,it);
            residuals(i_rows) = res(eq_index);
            if it == 2
                [rows,cols,vals] = find(jacobian(eq_index,i_cols_1));
                iJacobian{1} = [offset+rows, i_cols_J1(cols), vals];
            elseif it == T + 1
                [rows,cols,vals] = find(jacobian(eq_index,i_cols_T));
                iJacobian{T} = [offset+rows, i_cols_J(i_cols_T(cols)), vals];
            else
                [rows,cols,vals] = find(jacobian(eq_index,i_cols_j));
                iJacobian{it-1} = [offset+rows, i_cols_J(cols), vals];
                i_cols_J = i_cols_J + ny;
            end
            offset = offset + ny;
        end

        i_rows = i_rows + ny;
        i_cols = i_cols + ny;
    end

    if nargout == 2
        iJacobian = cat(1,iJacobian{:});
        JJacobian = sparse(iJacobian(:,1),iJacobian(:,2),iJacobian(:,3),T* ...
                           ny,T*ny);
    end