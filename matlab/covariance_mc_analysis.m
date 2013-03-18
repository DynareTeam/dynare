function oo_ = covariance_mc_analysis(NumberOfSimulations,type,dname,fname,vartan,nvar,var1,var2,mh_conf_sig,oo_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variables covariance matrix.

% Copyright (C) 2008-2009 Dynare Team
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

if strcmpi(type,'posterior')
    TYPE = 'Posterior';
    PATH = [dname '/metropolis/'];
else
    TYPE = 'Prior';
    PATH = [dname '/prior/moments/'];
end

indx1 = check_name(vartan,var1);
if isempty(indx1)
    disp([ type '_analysis:: ' var1 ' is not a stationary endogenous variable!'])
    return
end
if ~isempty(var2)
    indx2 = check_name(vartan,var2);
    if isempty(indx2)
        disp([ type '_analysis:: ' var2 ' is not a stationary endogenous variable!'])
        return
    end
else
    indx2 = indx1;
    var2 = var1;
end

if isfield(oo_,[ TYPE 'TheoreticalMoments'])
    eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments;'])   
    if isfield(temporary_structure,'dsge')
        eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge;'])
        if isfield(temporary_structure,'covariance')
            eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Mean;'])
            if isfield(temporary_structure,var1)
                eval(['temporary_structure_1 = oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Mean.' var1 ';'])
                if isfield(temporary_structure_1,var2)
                    % Nothing to do (the covariance matrix is symmetric!).
                    return
                end
            else
                if isfield(temporary_structure,var2)
                    eval(['temporary_structure_2 = oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Mean.' var2 ';'])
                    if isfield(temporary_structure_2,var1)
                        % Nothing to do (the covariance matrix is symmetric!).
                        return
                    end
                end
            end
        end
    end
end

ListOfFiles = dir([ PATH  fname '_' TYPE '2ndOrderMoments*.mat']);
i1 = 1; tmp = zeros(NumberOfSimulations,1);
for file = 1:length(ListOfFiles)
    load([ PATH ListOfFiles(file).name ]);
    i2 = i1 + rows(Covariance_matrix) - 1;
    tmp(i1:i2) = Covariance_matrix(:,symmetric_matrix_index(indx1,indx2,nvar));
    i1 = i2+1;
end
name = [var1 '.' var2];
if ~isconst(tmp)
    [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Mean.' name ' = p_mean;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Median.' name ' = p_median;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Variance.' name ' = p_var;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.HPDinf.' name ' = hpd_interval(1);']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.HPDsup.' name ' = hpd_interval(2);']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.deciles.' name ' = p_deciles;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.density.' name ' = density;']);
else
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Mean.' name ' = NaN;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Median.' name ' = NaN;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.Variance.' name ' = NaN;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.HPDinf.' name ' = NaN;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.HPDsup.' name ' = NaN;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.deciles.' name ' = NaN;']);
    eval(['oo_.' TYPE 'TheoreticalMoments.dsge.covariance.density.' name ' = NaN;']);
end