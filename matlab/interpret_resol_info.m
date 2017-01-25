function message = interpret_resol_info(info)
    
% Returns a message describing problem encountered during the resolution of
% a model.
%
% INPUTS 
% - info       [struct]  Second output argument return by the resol routine
%
% OUTPUTS 
% - message    [string]  Description of the issue preventing model's resolution.
 
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

switch info(1)
  case 0
    message = '';
  case 1
    message = 'The model doesn''t determine the current variable uniquely.';
  case 2
    message = 'MJDGGES (Generalized Schur decomposition) returned an error code.';
  case 3
    message = 'Blanchard & Kahn conditions are not satisfied: no stable equilibrium.';
  case 4
    message = 'Blanchard & Kahn conditions are not satisfied: indeterminacy.';
  case 5
    message = 'Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.';
  case 6
    message = 'The jacobian evaluated at the deterministic steady state is complex.';
  case 19
    message = 'The steadystate routine thrown an exception (inconsistent deep parameters).';
  case 20
    message = sprintf('Cannot find the steady state (the sum of square residuals of the static equations is %s)', num2str(info(2)));
  case 21
    message = sprintf('The steady state is complex (the sum of square residuals of imaginary parts of the steady state is %s)', num2str(info(2)));
  case 22
    message = 'The steady state has NaNs.';
  case 23
    message = 'Parameters have been updated in the steadystate routine and some have complex values.';
  case 24
    message = 'Parameters have been updated in the steadystate routine and some are NaNs.';
  case 30
    message = 'Ergodic variance can''t be computed.';
  otherwise
    message = 'Unknown issue!';
end