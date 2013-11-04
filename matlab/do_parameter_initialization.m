function [xparam1,estim_params_,xparam1_explicitly_initialized,xparam1_properly_calibrated]=do_parameter_initialization(estim_params_,xparam1_calib,xparam1_NaN_set_to_prior_mean)
% function [xparam1,estim_params_]=get_initialized_parameters(estim_params_,xparam1_calib)
% gets explicitly initialized variables and properly calibrated parameters
%
% INPUTS
%    o estim_params_    [structure] characterizing parameters to be estimated.
%    o xparam1_calib    [double]    vector of parameters to be estimated, with parameters
%                                   initialized from calibration using get_all_parameters
%     
%    o xparam1_NaN_set_to_prior_mean [double]    vector of parameters to be estimated, with parameters
%                                                initialized using dynare_estimation_init; not explicitly initialized
%                                                parameters are at prior mean
% OUTPUTS
%    o xparam1                           [double]    vector of initialized parameters; uses the hierarchy: 1) explicitly initialized parameters,
%                                                    2) calibrated parameters, 3) prior mean
%    o estim_params_    [structure] characterizing parameters to be estimated; it is 
%                                   updated here to reflect calibrated parameters
%    o xparam1_explicitly_initialized    [double]    vector of parameters to be estimated that
%                                                    were explicitly initialized 
%    o xparam1_properly_calibrated       [double]    vector of parameters to be estimated that
%                                                    were properly calibrated 
%    
% SPECIAL REQUIREMENTS
%    None

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

nvx = size(estim_params_.var_exo,1);
nvn = size(estim_params_.var_endo,1);
ncx = size(estim_params_.corrx,1);
ncn = size(estim_params_.corrn,1);
np = size(estim_params_.param_vals,1);

estim_params_.nvx = nvx; %exogenous shock variances
estim_params_.nvn = nvn; %endogenous variances, i.e. measurement error
estim_params_.ncx = ncx; %exogenous shock correlations
estim_params_.ncn = ncn; % correlation between endogenous variables, i.e. measurement error.
estim_params_.np = np;   % other parameters of the model

xparam1_explicitly_initialized = NaN(nvx+nvn+ncx+ncn+np,1);
xparam1_properly_calibrated = NaN(nvx+nvn+ncx+ncn+np,1);

offset=0;
if nvx
    initialized_par_index=find(~isnan(estim_params_.var_exo(:,2)));
    calibrated_par_index=find(isnan(estim_params_.var_exo(:,2)) & ~isnan(xparam1_calib(offset+1:offset+nvx,1)));
    uninitialized_par_index=find(isnan(estim_params_.var_exo(:,2)) & isnan(xparam1_calib(offset+1:offset+nvx,1)));    
    xparam1_explicitly_initialized(offset+initialized_par_index,1) = estim_params_.var_exo(initialized_par_index,2);
    %update estim_params_ with calibrated starting values
    estim_params_.var_exo(calibrated_par_index,2)=xparam1_calib(offset+calibrated_par_index,1);
    %find parameters that are calibrated and do not violate inverse gamma prior
    xparam1_properly_calibrated(offset+calibrated_par_index,1) = xparam1_calib(offset+calibrated_par_index,1);
    inv_gamma_violation=find(estim_params_.var_exo(calibrated_par_index,2)==0 & estim_params_.var_exo(calibrated_par_index,4)==0);
    if inv_gamma_violation
        estim_params_.var_exo(calibrated_par_index(inv_gamma_violation),2)=NaN;        
        xparam1_properly_calibrated(offset+calibrated_par_index(inv_gamma_violation),1)=NaN;        
        fprintf('PARAMETER INITIALIZATION: Some standard deviations of shocks of the calibrated model are 0 and\n') 
        fprintf('PARAMETER INITIALIZATION: violate the inverse gamma prior. They will instead be initialized with the prior mean.\n')    
    end
    if uninitialized_par_index
        fprintf('PARAMETER INITIALIZATION: Warning, some estimated standard deviations of shocks are not\n') 
        fprintf('PARAMETER INITIALIZATION: initialized. They will be initialized with the prior mean.\n')
    end
end
offset=offset+nvx;
if nvn
    initialized_par_index=find(~isnan(estim_params_.var_endo(:,2)));
    calibrated_par_index=find(isnan(estim_params_.var_endo(:,2)) & ~isnan(xparam1_calib(offset+1:offset+nvn,1)));
    uninitialized_par_index=find(isnan(estim_params_.var_endo(:,2)) & isnan(xparam1_calib(offset+1:offset+nvn,1)));
    xparam1_explicitly_initialized(offset+initialized_par_index,1) = estim_params_.var_endo(initialized_par_index,2);
    estim_params_.var_endo(calibrated_par_index,2)=xparam1_calib(offset+calibrated_par_index,1);
    %find parameters that are calibrated and do not violate inverse gamma prior
    xparam1_properly_calibrated(offset+calibrated_par_index,1) = xparam1_calib(offset+calibrated_par_index,1);
    inv_gamma_violation=find(estim_params_.var_endo(calibrated_par_index,2)==0 & estim_params_.var_endo(calibrated_par_index,5)==4);
    if inv_gamma_violation
        estim_params_.var_endo(calibrated_par_index(inv_gamma_violation),2)=NaN;        
        xparam1_properly_calibrated(offset+calibrated_par_index(inv_gamma_violation),1)=NaN;        
        fprintf('PARAMETER INITIALIZATION: Some measurement errors of the calibrated model are 0 and violate the\n')
        fprintf('PARAMETER INITIALIZATION: inverse gamma prior. They will instead be initialized with the prior mean.\n')    
    end
    if uninitialized_par_index
        fprintf('PARAMETER INITIALIZATION: Warning, some measurement errors are not initialized. They will be initialized\n')
        fprintf('PARAMETER INITIALIZATION: with the prior mean.\n')
    end
end
offset=offset+nvn;
if ncx
    initialized_par_index=find(~isnan(estim_params_.corrx(:,3)));
    calibrated_par_index=find(isnan(estim_params_.corrx(:,3)) & ~isnan(xparam1_calib(offset+1:offset+ncx,1)));
    uninitialized_par_index=find(isnan(estim_params_.corrx(:,3)) & isnan(xparam1_calib(offset+1:offset+ncx,1)));
    xparam1_explicitly_initialized(offset+initialized_par_index,1) = estim_params_.corrx(initialized_par_index,3);
    estim_params_.corrx(calibrated_par_index,3)=xparam1_calib(offset+calibrated_par_index,1);
    xparam1_properly_calibrated(offset+calibrated_par_index,1) = xparam1_calib(offset+calibrated_par_index,1);
    if uninitialized_par_index
        fprintf('PARAMETER INITIALIZATION: Warning, some correlations between structural shocks are not initialized.\n')
        fprintf('PARAMETER INITIALIZATION: They will be initialized with the prior mean.\n')
    end
end
offset=offset+ncx;
if ncn
    initialized_par_index=find(~isnan(estim_params_.corrn(:,3)));
    calibrated_par_index=find(isnan(estim_params_.corrn(:,3)) & ~isnan(xparam1_calib(offset+1:offset+ncn,1)));
    uninitialized_par_index=find(isnan(estim_params_.corrn(:,3)) & isnan(xparam1_calib(offset+1:offset+ncn,1)));
    xparam1_explicitly_initialized(offset+initialized_par_index,1) = estim_params_.corrn(initialized_par_index,3);
    estim_params_.corrn(calibrated_par_index,3)=xparam1_calib(offset+calibrated_par_index,1);
    xparam1_properly_calibrated(offset+calibrated_par_index,1) = xparam1_calib(offset+calibrated_par_index,1);
    if uninitialized_par_index
        fprintf('PARAMETER INITIALIZATION: Warning, some correlations between measurement errors are not initialized.\n')
        fprintf('PARAMETER INITIALIZATION: They will be initialized with the prior mean.\n')
    end
end
offset=offset+ncn;
if np
    initialized_par_index=find(~isnan(estim_params_.param_vals(:,2)));
    calibrated_par_index=find(isnan(estim_params_.param_vals(:,2)) & ~isnan(xparam1_calib(offset+1:offset+np,1)));
    uninitialized_par_index=find(isnan(estim_params_.param_vals(:,2)) & isnan(xparam1_calib(offset+1:offset+np,1)));
    xparam1_explicitly_initialized(offset+initialized_par_index,1) = estim_params_.param_vals(initialized_par_index,2);
    estim_params_.param_vals(calibrated_par_index,2)=xparam1_calib(offset+calibrated_par_index,1);
    xparam1_properly_calibrated(offset+calibrated_par_index,1) = xparam1_calib(offset+calibrated_par_index,1);
    if uninitialized_par_index
        fprintf('PARAMETER INITIALIZATION: Warning, some deep parameters are not initialized. They will be\n')
        fprintf('PARAMETER INITIALIZATION: initialized with the prior mean.\n')
    end
end
xparam1=xparam1_explicitly_initialized;
xparam1(isnan(xparam1))=xparam1_properly_calibrated(isnan(xparam1)); %set not explicitly initialized parameters that do not obviously violate prior distribution to calibrated parameter values
xparam1(isnan(xparam1))=xparam1_NaN_set_to_prior_mean(isnan(xparam1)); %set not yet initialized parameters to prior mean coming from dynare_estimation_init

