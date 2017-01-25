function [fval,exit_flag,arg1,arg2] = penalty_objective_function(x0,fcn,penalty,varargin)
%function [fval,exit_flag,arg1,arg2] = penalty_objective_function(x0,fcn,penalty,varargin)

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

    [fval,info,exit_flag,arg1,arg2] = fcn(x0,varargin{:});

    if info(1) ~= 0
        fval = penalty + info(4);
    end
end