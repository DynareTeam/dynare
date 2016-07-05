function [res,fjac,domer] = mcp_func(x,jacflag)
% function [res,fjac,domer] = mcp_func(x,jacflag)
% wrapper function for mixed complementarity problem when using PATH
% 
% INPUTS 
% - x                   [double] N*T array, paths for the endogenous variables (initial guess).
% - jacflag             [scalar] indicator whether Jacobian is requested
%
% OUTPUTS 
%  - res                [double] (N*T)*1 array, residuals of the stacked problem
%  - fjac               [double] (N*T)*(N*T) array, Jacobian of the stacked problem
%  - domer              [scalar] errorflag that is 1 if solution is not real

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

global mcp_data

if jacflag
    [res,fjac] = mcp_data.func(x,mcp_data.args{:});
    fjac = sparse(fjac);
else
    res = mcp_data.func(x,mcp_data.args{:});
    fjac = [];
end
if isreal(res)
    domer = 0;
else
    domer = 1;
end
