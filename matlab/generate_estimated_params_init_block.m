function generate_estimated_params_init_block(xparam1,estim_params_,M_,options_)

% function generate_estimated_params_init_block(xparam1,estim_params_,M_,options_)
% writes parameter values from xparam_calib derived from get_all_parameters into 
% an estimated_params_init-block that can be used to start the estimation.
% 
% INPUTS
%    xparam1:        Parameter vector from which to create block
%	 estim_params_:  Dynare structure describing the estimated parameters.
%    M_:             Dynare structure describing the model. 
%    options_:       Dynare options structure
%
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2013 Dynare Team
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

nvx = estim_params_.nvx;
ncx = estim_params_.ncx;
nvn = estim_params_.nvn;
ncn = estim_params_.ncn;
np = estim_params_.np;
% stderrs of the exogenous shocks
fprintf('\nestimated_params_init;\n')
if nvx
    for ii=1:nvx
        vname = deblank(M_.exo_names(estim_params_.var_exo(ii,1),:));
        fprintf('stderr %s, %f;\n', vname,xparam1(ii));
    end
end
% update offset
offset = nvx;

% setting measument error variance
if nvn
    for ii=1:nvn
        vname = deblank(options_.varobs(estim_params_.nvn_observable_correspondence(ii,1),:));
        fprintf('stderr %s, %f;\n', vname,xparam1(offset+ii));
    end
end

% update offset
offset = nvx+nvn;

% correlations among shocks (ncx)
if ncx
    corrx = estim_params_.corrx;
    for ii=1:ncx
        k1 = corrx(ii,1);
        k2 = corrx(ii,2);
        vname1 = deblank(M_.exo_names(k1,:)); 
        vname2 = deblank(M_.exo_names(k2,:));
        fprintf('corr %s, %s, %f;\n', vname1,vname2,xparam1(offset+ii));
    end
end
% update offset
offset = nvx+nvn+ncx;

if ncn
    for ii=1:ncn
        vname1 = deblank(options_.varobs(estim_params_.corrn_observable_correspondence(ii,1),:));
        vname2 = deblank(options_.varobs(estim_params_.corrn_observable_correspondence(ii,2),:));
        fprintf('corr %s, %s, %f;\n', vname1,vname2,xparam1(offset+ii));
    end
end

% update offset
offset = nvx+ncx+nvn+ncn;


% structural parameters
if np
    for ii=1:np
        jj1 = estim_params_.param_vals(ii,1);
        vname = deblank(M_.param_names(jj1,:));
        fprintf('%s, %f;\n',vname,xparam1(offset+ii));
    end 
end
fprintf('end;\n')
