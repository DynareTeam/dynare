function info=test_for_deep_parameters_calibration(M_)
% Issues a warning is some of the parameters are NaNs.
%
% INPUTS
%   M_    [structure]   Description of the (simulated or estimated) model.
%  
% OUTPUTS
%   info  [scalar]      0 if no problems detected, 1 otherwise 
%    
% ALGORITHM
%   none
%    
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2010-2014 Dynare Team
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
plist = list_of_parameters_calibrated_as_NaN(M_);
if ~isempty(plist)
    info=1;
    message = ['Some of the parameters have no value (' ];
    for i=1:size(plist,1)
        if i<size(plist,1)
            message = [message, deblank(plist(i,:)) ', '];
        else
            message = [message, deblank(plist(i,:)) ')'];
        end
    end
    tmp = dbstack;
    message = [message, ' when using ' tmp(2).name '. '];
    message = [message, 'If these parameters are not initialized in a steadystate file or a steady_state_model-block, Dynare may not be able to solve the model...'];
    message_id  = 'Dynare:ParameterCalibration:NaNValues';
    warning(message_id,message);
    if strmatch('optimal_policy_discount_factor',plist,'exact')
        warning('Either you have not correctly initialized planner_discount or you are calling a command like steady or stoch_simul that is not allowed in the context of ramsey_policy')
    end
else
    info=0;
end