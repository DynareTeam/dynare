function [steady_state,params,info] = steady_(M_,options_,oo_)
% function [steady_state,params,info] = steady_(M_,options_,oo_)
% Computes the steady state 
%  
% INPUTS
%   M                         struct           model structure
%   options                   struct           options
%   oo                        struct           output results
%  
% OUTPUTS
%   steady_state              vector           steady state
%   params                    vector           parameters (may have been
%                                              modified by user in
%                                              explicit computation of
%                                              the steady state)
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2016 Dynare Team
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

if options_.solve_algo < 0 || options_.solve_algo > 12
    error('STEADY: solve_algo must be between 0 and 12')
end

if ~options_.bytecode && ~options_.block && options_.solve_algo > 4 && ...
        options_.solve_algo < 10
    error('STEADY: you can''t use solve_algo > 4 without block nor bytecode options')
end

if ~options_.bytecode && options_.block && options_.solve_algo == 5
    error('STEADY: you can''t use solve_algo = 5 without bytecode option')
end

if isoctave && ismember(options_.solve_algo,[7,11])
    error(['SIMUL: you can''t use solve_algo = %u under Octave'],options_.solve_algo)
end

[steady_state,params,info] = evaluate_steady_state(oo_.steady_state,M_,options_,oo_,~options_.steadystate.nocheck);