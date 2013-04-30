function i = name2index(options_, M_, estim_params_, type, name1, name2 )
% Returns the index associated to an estimated object (deep parameter,
% variance of a structural shock or measurement error, covariance between
% two structural shocks, covariance between two measurement errors).
%  
% INPUTS:
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   estim_params_   [structure]    Dynare structure (related to estimation).
%   type            [string]       'DeepParameter', 'MeasurementError' (for measurement equation error) or 'StructuralShock' (for structural shock).
%   name1           [string]       
%   name2           [string]    
% OUTPUTS
%   i               [integer]      column index (in x2, an array of mh draws) associated to name1[,name2].
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2008-2013 Dynare Team
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

nvx     = estim_params_.nvx;
nvn     = estim_params_.nvn;
ncx     = estim_params_.ncx;
ncn     = estim_params_.ncn;
npa     = estim_params_.np ;
nnn = nvx+nvn+ncx+ncn+npa;

i = [];

if strcmpi(type,'DeepParameter')
    i = nvx + nvn + ncx + ncn + ...
        strmatch(name1,M_.param_names(estim_params_.param_vals(:,1),:),'exact');
    if nargin>5
        disp('The last input argument is useless!')
    end
    if isempty(i)
        disp([name1 ' is not an estimated deep parameter!'])
    end
    return
end

if strcmpi(type,'StructuralShock')
    if nargin<6% Covariance matrix diagonal term.
        i = strmatch(name1,M_.exo_names(estim_params_.var_exo(:,1),:),'exact');
        if isempty(i)
            disp(['The standard deviation of ' name1  ' is not an estimated parameter!'])
            return
        end
    else% Covariance matrix off-diagonal term
        offset = nvx+nvn;
        try 
            list_of_structural_shocks = { M_.exo_names(estim_params_.corrx(:,1),:) , M_.exo_names(estim_params_.corrx(:,2),:) };
            k1 = strmatch(name1,list_of_structural_shocks(:,1),'exact');
            k2 = strmatch(name2,list_of_structural_shocks(:,2),'exact');
            i = offset+intersect(k1,k2);
            if isempty(i)
                k1 = strmatch(name1,list_of_structural_shocks(:,2),'exact');
                k2 = strmatch(name2,list_of_structural_shocks(:,1),'exact');
                i = offset+intersect(k1,k2);
            end
            if isempty(i)
                if isempty(i)
                    disp(['The correlation between ' name1 ' and ' name2 ' is not an estimated parameter!'])
                    return
                end
            end
        catch
            disp(['Off diagonal terms of the covariance matrix are not estimated (state equation)'])
        end
    end
end

if strcmpi(type,'MeasurementError')
    if nargin<6% Covariance matrix diagonal term
        i = nvx + strmatch(name1,M_.endo_names(estim_params_.var_endo(:,1),:),'exact');
        if isempty(i)
            disp(['The standard deviation of the measurement error on ' name1  ' is not an estimated parameter!'])
            return
        end
    else% Covariance matrix off-diagonal term
        offset = nvx+nvn+ncx;
        try
            list_of_measurement_errors = { M_.endo_names(estim_params_.corrn(:,1),:) , M_.endo_names(estim_params_.corrn(:,2),:) };
            k1 = strmatch(name1,list_of_measurement_errors(:,1),'exact');
            k2 = strmatch(name2,list_of_measurement_errors(:,2),'exact');
            i = offset+intersect(k1,k2);
            if isempty(i)
                k1 = strmatch(name1,list_of_measurement_errors(:,2),'exact');
                k2 = strmatch(name2,list_of_measurement_errors(:,1),'exact');
                i = offset+intersect(k1,k2);
            end
            if isempty(i)
                if isempty(i)
                    disp(['The correlation between the measurement errors on ' name1 ' and ' name2 ' is not an estimated parameter!'])
                    return
                end
            end
        catch
            disp('Off diagonal terms of the covariance matrix are not estimated (measurement equation)')
        end
    end
end