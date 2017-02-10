function [fval, exit_flag, arg1, arg2] = penalty_objective_function(x, fcn, base_penalty, varargin)

% Encapsulates an objective function to be minimized, adding a penalty if necessary.
%
% INPUTS 
% - x             [double]    n*1 vector of instrument values.
% - fcn           [fhandle]   objective function.
% - base_penalty  [double]    scalar, base of the penality (typically the value of the objective at the previous iteration).
% - varagin       [cell]      additional parameters for fcn.
%
% OUTPUTS 
% - fval          [double]    scalar, value of the objective function at x.
% - exit_flag     [integer]   scalar, flag returned by fcn (third output).
% - arg1, arg2                fourth and fifth output arguments of the objective function.

% Copyright (C) 2016-2017 Dynare Team
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

[fval, info, exit_flag, arg1, arg2] = fcn(x, varargin{:});

if info(1)
    fval = base_penalty + info(4);
end