function xparam = get_posterior_parameters(type,M_,estim_params_,oo_,options_,field1)

% function xparam = get_posterior_parameters(type,M_,estim_params_,oo_,options_,field1)
% Selects (estimated) parameters (posterior mode or posterior mean).
%
% INPUTS
%   o type              [char]     = 'mode' or 'mean'.
%   o M_:               [structure] Dynare structure describing the model.
%   o estim_params_:    [structure] Dynare structure describing the estimated parameters.
%   o field_1           [char]     optional field like 'mle_'.
%
% OUTPUTS
%   o xparam     vector of estimated parameters
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2017 Dynare Team
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

if nargin<6
    field1='posterior_';
end
nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np;

xparam = zeros(nvx+nvn+ncx+ncn+np,1);

m = 1;
for i=1:nvx
    k1 = estim_params_.var_exo(i,1);
    name1 = deblank(M_.exo_names(k1,:));
    xparam(m) = oo_.([field1 type]).shocks_std.(name1);
    m = m+1;
end

for i=1:nvn
    k1 = estim_params_.nvn_observable_correspondence(i,1);
    name1 = options_.varobs{k1};
    xparam(m) = oo_.([field1 type]).measurement_errors_std.(name1);
    m = m+1;
end

for i=1:ncx
    k1 = estim_params_.corrx(i,1);
    k2 = estim_params_.corrx(i,2);
    name1 = deblank(M_.exo_names(k1,:));
    name2 = deblank(M_.exo_names(k2,:));
    xparam(m) = oo_.([field1 type]).shocks_corr.([name1 '_' name2]);
    m = m+1;
end

for i=1:ncn
    k1 = estim_params_.corrn_observable_correspondence(i,1);
    k2 = estim_params_.corrn_observable_correspondence(i,2);
    name1 = options_.varobs{k1};
    name2 = options_.varobs{k2};
    xparam(m) = oo_.([field1 type]).measurement_errors_corr.([name1 '_' name2]);
    m = m+1;
end

FirstDeep = m;

for i=1:np
    name1 = deblank(M_.param_names(estim_params_.param_vals(i,1),:));
    xparam(m) = oo_.([field1 type]).parameters.(name1);
    m = m+1;
end