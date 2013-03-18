function oo_ = variance_decomposition_mc_analysis(NumberOfSimulations,type,dname,fname,exonames,exo,vartan,var,mh_conf_sig,oo_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variance decomposition.

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

indx = check_name(vartan,var);
if isempty(indx)
    disp([ type '_analysis:: ' var ' is not a stationary endogenous variable!'])
    return
end
jndx = check_name(exonames,exo);
if isempty(jndx)
    disp([ type '_analysis:: ' exo ' is not a declared exogenous variable!'])
    return
end

name = [ var '.' exo ];
if isfield(oo_, [ TYPE 'TheoreticalMoments'])
    eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments;'])
    if isfield(temporary_structure,'dsge')
        eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge;'])
        if isfield(temporary_structure,'VarianceDecomposition')
            eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.Mean;'])
            if isfield(temporary_structure,name)
                % Nothing to do.
                return
            end
        end
    end
end

ListOfFiles = dir([ PATH  fname '_' TYPE 'VarianceDecomposition*.mat']);
i1 = 1; tmp = zeros(NumberOfSimulations,1);
indice = (indx-1)*rows(exonames)+jndx;
for file = 1:length(ListOfFiles)
    load([ PATH ListOfFiles(file).name ]);
    i2 = i1 + rows(Decomposition_array) - 1;
    tmp(i1:i2) = Decomposition_array(:,indice);
    i1 = i2+1;
end

t1 = min(tmp); t2 = max(tmp);
t3 = t2-t1;% How to normalize ? t1 and t2 may be zero...
if t3<1.0e-12
    if t1<1.0e-12
        t1 = 0;
    end
    if abs(t1-1)<1.0e-12
        t1 = 1;
    end 
    p_mean = t1;
    p_median = t1;
    p_var = 0;
    hpd_interval = NaN(2,1);
    p_deciles = NaN(9,1);
    density = NaN;
else
    [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
end
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.Mean.' name ' = p_mean;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.Median.' name ' = p_median;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.Variance.' name ' = p_var;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.HPDinf.' name ' = hpd_interval(1);']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.HPDsup.' name ' = hpd_interval(2);']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.deciles.' name ' = p_deciles;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.VarianceDecomposition.density.' name ' = density;']);