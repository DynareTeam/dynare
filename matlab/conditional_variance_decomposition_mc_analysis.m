function oo_ = ...
    conditional_variance_decomposition_mc_analysis(NumberOfSimulations, type, dname, fname, Steps, exonames, exo, var_list, endogenous_variable_index, mh_conf_sig, oo_)
% This function analyses the (posterior or prior) distribution of the
% endogenous conditional variance decomposition.

% Copyright (C) 2009-2010 Dynare Team
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

% $$$ indx = check_name(vartan,var);
% $$$ if isempty(indx)
% $$$     disp([ type '_analysis:: ' var ' is not a stationary endogenous variable!'])
% $$$     return
% $$$ end
% $$$ endogenous_variable_index = sum(1:indx);
exogenous_variable_index = check_name(exonames,exo);
if isempty(exogenous_variable_index)
    disp([ type '_analysis:: ' exo ' is not a declared exogenous variable!'])
    return
end

name = [ var_list(endogenous_variable_index,:) '.' exo ];
if isfield(oo_, [ TYPE 'TheoreticalMoments' ])
    eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments;'])
    if isfield(temporary_structure,'dsge')
        eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge;'])
        if isfield(temporary_structure,'ConditionalVarianceDecomposition')
            eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean;'])
            if isfield(temporary_structure,name)
                if sum(Steps-temporary_structure.(name)(1,:)) == 0
                    % Nothing (new) to do here...
                    return
                end
            end
        end
    end
end

ListOfFiles = dir([ PATH  fname '_' TYPE 'ConditionalVarianceDecomposition*.mat']);
i1 = 1; tmp = zeros(NumberOfSimulations,length(Steps));
for file = 1:length(ListOfFiles)
    load([ PATH ListOfFiles(file).name ]);
    % 4D-array (endovar,time,exovar,simul)
    i2 = i1 + size(Conditional_decomposition_array,4) - 1;
    tmp(i1:i2,:) = transpose(dynare_squeeze(Conditional_decomposition_array(endogenous_variable_index,:,exogenous_variable_index,:)));
    i1 = i2+1;
end

p_mean = NaN(1,length(Steps));
p_median = NaN(1,length(Steps));
p_variance = NaN(1,length(Steps));
p_deciles = NaN(9,length(Steps));
p_density = NaN(2^9,2,length(Steps));
p_hpdinf = NaN(1,length(Steps));
p_hpdsup = NaN(1,length(Steps));
for i=1:length(Steps)
    [pp_mean, pp_median, pp_var, hpd_interval, pp_deciles, pp_density] = ...
        posterior_moments(tmp(:,i),1,mh_conf_sig);
    p_mean(i) = pp_mean;
    p_median(i) = pp_median;
    p_variance(i) = pp_var;
    p_deciles(:,i) = pp_deciles;
    p_hpdinf(i) = hpd_interval(1);
    p_hpdsup(i) = hpd_interval(2);
    p_density(:,:,i) = pp_density;
end
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.Steps = Steps;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean.' name ' = p_mean;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.Median.' name ' = p_median;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.Variance.' name ' = p_variance;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.HPDinf.' name ' = p_hpdinf;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.HPDsup.' name ' = p_hpdsup;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.deciles.' name ' = p_deciles;']);
eval(['oo_.' TYPE 'TheoreticalMoments.dsge.ConditionalVarianceDecomposition.density.' name ' = p_density;']);