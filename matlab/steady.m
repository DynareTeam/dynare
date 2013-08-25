function steady()
% function steady()
% computes and prints the steady state calculations
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2012 Dynare Team
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

global M_ oo_ options_ ys0_ 

test_for_deep_parameters_calibration(M_);

if options_.steadystate_flag && options_.homotopy_mode
    error('STEADY: Can''t use homotopy when providing a steady state external file');
end


% Keep of a copy of M_.Sigma_e
Sigma_e = M_.Sigma_e;

% Set M_.Sigma_e=0 (we compute the *deterministic* steady state)
M_.Sigma_e = zeros(size(Sigma_e));

info = 0;
switch options_.homotopy_mode
  case 1
    [M_,oo_,info,ip,ix,ixd] = homotopy1(options_.homotopy_values,options_.homotopy_steps,M_,options_,oo_);
  case 2
    homotopy2(options_.homotopy_values, options_.homotopy_steps);
  case 3
    [M_,oo_,info,ip,ix,ixd] = homotopy3(options_.homotopy_values,options_.homotopy_steps,M_,options_,oo_);
end

if info(1)
    hv = options_.homotopy_values;
    skipline()
    disp('WARNING: homotopy step was not completed')
    disp('The last values for which a solution was found are:')
    for i=1:length(ip)
        disp(sprintf('%12s %12.6f',M_.param_names(hv(ip(i),2),:), ...
                     M_.params(hv(ip(i),2))))
    end
    for i=1:length(ix)
        disp(sprintf('%12s %12.6f',M_.exo_names(hv(ix(i),2),:), ...
                     oo_.exo_steady_state(hv(ix(i),2))))
    end
    for i=1:length(ixd)
        disp(sprintf('%12s %12.6f',M_.exo_det_names(hv(ixd(i),2),:), ...
                     oo_.exo_det_steady_state(hv(ixd(i),2))))
    end
    
    if options_.homotopy_force_continue
        disp('Option homotopy_continue is set, so I continue ...')
    else
        error('Homotopy step failed')
    end
end
    
[steady_state,M_.params,info] = steady_(M_,options_,oo_);
oo_.steady_state = steady_state;

if info(1) == 0
    if options_.noprint == 0
        disp_steady_state(M_,oo_);
    end
else
    if options_.noprint == 0
        resid;
    end
    print_info(info,options_.noprint, options_);
end

M_.Sigma_e = Sigma_e;
