function xparam1=get_all_parameters(estim_params_,M_)

% function xparam1=get_parameters
% gets parameters values from M_.params into xparam1 (inverse mapping to set_all_parameters)
% This is called if a model was calibrated before estimation to back out
% parameter values
% 
% INPUTS
%	 estim_params_:  Dynare structure describing the estimated parameters.
%    M_:             Dynare structure describing the model. 
%
% OUTPUTS
%    xparam1:       N*1 double vector of parameters from calibrated model that are to be estimated
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
Sigma_e = M_.Sigma_e;
Correlation_matrix = M_.Correlation_matrix;
H = M_.H;
Correlation_matrix_ME = M_.Correlation_matrix_ME;

xparam1=NaN(nvx+ncx+nvn+ncn+np,1);
% stderrs of the exogenous shocks
if nvx
    var_exo = estim_params_.var_exo;
    for i=1:nvx
        k = var_exo(i,1);
        xparam1(i)=sqrt(Sigma_e(k,k));
    end
end
% update offset
offset = nvx;

% setting measument error variance
if nvn
    for i=1:nvn
        k = estim_params_.nvn_observable_correspondence(i,1);
        xparam1(offset+i)=sqrt(H(k,k));
    end
end

% update offset
offset = nvx+nvn;

% correlations among shocks (ncx)
if ncx
    corrx = estim_params_.corrx;
    for i=1:ncx
        k1 = corrx(i,1);
        k2 = corrx(i,2);
        xparam1(i+offset)=Correlation_matrix(k1,k2);
    end
end
% update offset
offset = nvx+nvn+ncx;

if ncn
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        xparam1(i+offset)=Correlation_matrix_ME(k1,k2);
    end
end

% update offset
offset = nvx+ncx+nvn+ncn;


% structural parameters
if np
     xparam1(offset+1:end)=M_.params(estim_params_.param_vals(:,1));
end