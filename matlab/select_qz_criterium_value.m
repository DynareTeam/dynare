function options_=select_qz_criterium_value(options_)
% function options_=select_qz_criterium_value(options_)
% set the value of options_.qz_criterium depending on the Kalman filter used
%
% INPUTS
%   options_:   Dynare options structure
%
% OUTPUTS
%   options_:   Dynare options structure
%
% SPECIAL REQUIREMENTS
%   none

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


% If options_.lik_init == 1
%     set by default options_.qz_criterium to 1-1e-6
%     and check options_.qz_criterium < 1-eps if options_.lik_init == 1
% Else
%     set by default options_.qz_criterium to 1+1e-6
stack = dbstack;

if isequal(options_.lik_init,1)
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1-1e-6;
    elseif options_.qz_criterium > 1-eps
        error([stack(2).file ': option qz_criterium is too large for estimating/smoothing ' ...
               'a stationary model. If your model contains unit roots, use ' ...
               'option diffuse_filter'])
    end
else
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    else
        if options_.qz_criterium <= 1
            fprintf('\n%s:: diffuse filter is incompatible with a qz_criterium<=1. Resetting it to 1+1e-6.\n',stack(2).file)
            options_.qz_criterium = 1+1e-6;
        end
    end
end
