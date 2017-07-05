function M_=set_parameters_locally(M_,xparam1)

% function M_out=set_parameters(M_,xparam1)
% Sets parameters value (except measurement errors)
% This is called for computations such as IRF and forecast
% when measurement errors aren't taken into account; in contrast to 
% set_parameters.m, the global M_-structure is not altered
%
% INPUTS
%    xparam1:   vector of parameters to be estimated (initial values)
%    M_:        Dynare model-structure
%
% OUTPUTS
%    M_:        Dynare model-structure
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017 Dynare Team
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

global estim_params_

nvx = estim_params_.nvx;
ncx = estim_params_.ncx;
nvn = estim_params_.nvn;
ncn = estim_params_.ncn;
np = estim_params_.np;
Sigma_e = M_.Sigma_e;
Correlation_matrix = M_.Correlation_matrix;
offset = 0;

% setting shocks variance on the diagonal of Covariance matrix; used later
% for updating covariances
if nvx
    var_exo = estim_params_.var_exo;
    for i=1:nvx
        k = var_exo(i,1);
        Sigma_e(k,k) = xparam1(i)^2;
    end
end
% and update offset
offset = offset + nvx + nvn;

% correlations amonx shocks (ncx)
if ncx
    corrx = estim_params_.corrx;
    for i=1:ncx
        k1 = corrx(i,1);
        k2 = corrx(i,2);
        Correlation_matrix(k1,k2) = xparam1(i+offset);
        Correlation_matrix(k2,k1) = Correlation_matrix(k1,k2);
    end
end
%build covariance matrix from correlation matrix and variances already on
%diagonal
Sigma_e = diag(sqrt(diag(Sigma_e)))*Correlation_matrix*diag(sqrt(diag(Sigma_e)));
if isfield(estim_params_,'calibrated_covariances')
    Sigma_e(estim_params_.calibrated_covariances.position)=estim_params_.calibrated_covariances.cov_value;
end

% and update offset
offset = offset + ncx + ncn;

% structural parameters
if np
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
end

M_.Sigma_e = Sigma_e;
M_.Correlation_matrix=Correlation_matrix;