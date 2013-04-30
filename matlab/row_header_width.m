function w=row_header_width(M_,estim_params_,bayestopt_)
% This function computes the width of the row headers for
% the estimation results
% 
% INPUTS 
%   estim_params_    [structure] 
%   M_               [structure]
%   bayestopt_       [structure]
%  
% OUTPUTS 
%   w                integer
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2013 Dynare Team
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

np = estim_params_.np;
nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;

w = 0;
if np
    for i=1:np
        w = max(w,length(bayestopt_.name{i}));
    end
end
if nvx
    for i=1:nvx
        k = estim_params_.var_exo(i,1); 
        w = max(w,length(deblank(M_.exo_names(k,:))));
    end
end
if nvn
    for i=1:nvn
        k = estim_params_.var_endo(i,1); 
        w = max(w,length(deblank(M_.endo_names(k,:))));
    end
end
if ncx
    for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        w = max(w,length(deblank(M_.exo_names(k1,:)))...
                +length(deblank(M_.exo_names(k2,:))));

    end
end
if ncn
    for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        w = max(w,length(deblank(M_.endo_names(k1,:)))...
                +length(deblank(M_.endo_names(k2,:))));

    end
end

