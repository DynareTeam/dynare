function oo_ = variance_decomposition_ME_mc_analysis(NumberOfSimulations,type,dname,fname,exonames,exo,vartan,var,mh_conf_sig,oo_,options_)
% function oo_ = variance_decomposition_ME_mc_analysis(NumberOfSimulations,type,dname,fname,exonames,exo,vartan,var,mh_conf_sig,oo_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variables' variance decomposition.
%
% INPUTS
%   NumberOfSimulations     [integer]           scalar, number of simulations.
%   type                    [string]            'prior' or 'posterior'
%   dname                   [string]            directory name where to save
%   fname                   [string]            name of the mod-file
%   exonames                [string]            (n_exo*char_length) character array with names of exogenous variables
%   exo                     [string]            name of current exogenous
%                                               variable
%   vartan                  [string]            (n_endo*char_length) character array with name
%                                               of endogenous variables
%   var                     [integer]         index of the current
%                                               endogenous variable
%   mh_conf_sig             [double]            2 by 1 vector with upper
%                                               and lower bound of HPD intervals
%   oo_                     [structure]         Dynare structure where the results are saved.
%   options_                [structure]         Dynare options structure
%
% OUTPUTS
%   oo_          [structure]        Dynare structure where the results are saved.



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
    if isequal(exo,'ME')
       jndx=size(exonames,1)+1; 
    else
        disp([ type '_analysis:: ' exo ' is not a declared exogenous variable!'])
    return
    end
end

var=deblank(var);
exo=deblank(exo);

name = [ var '.' exo ];
if isfield(oo_, [ TYPE 'TheoreticalMoments'])
    temporary_structure = oo_.([TYPE, 'TheoreticalMoments']);
    if isfield(temporary_structure,'dsge')
        temporary_structure = oo_.([TYPE, 'TheoreticalMoments']).dsge;
        if isfield(temporary_structure,'VarianceDecompositionME')
            temporary_structure = oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.Mean;
            if isfield(temporary_structure,name)
                % Nothing to do.
                return
            end
        end
    end
end

ListOfFiles = dir([ PATH  fname '_' TYPE 'VarianceDecompME*.mat']);
i1 = 1; tmp = zeros(NumberOfSimulations,1);
indice = (indx-1)*rows(exonames)+jndx;
for file = 1:length(ListOfFiles)
    load([ PATH ListOfFiles(file).name ]);
    i2 = i1 + rows(Decomposition_array_ME) - 1;
    tmp(i1:i2) = Decomposition_array_ME(:,indice);
    i1 = i2+1;
end

if options_.estimation.moments_posterior_density.indicator
    [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
else
    [p_mean, p_median, p_var, hpd_interval, p_deciles] = ...
        posterior_moments(tmp,0,mh_conf_sig);
end

oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.Mean.(var).(exo) = p_mean;
oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.Median.(var).(exo) = p_median;
oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.Variance.(var).(exo) = p_var;
oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.HPDinf.(var).(exo) = hpd_interval(1);
oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.HPDsup.(var).(exo) = hpd_interval(2);
oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.deciles.(var).(exo) = p_deciles;
if options_.estimation.moments_posterior_density.indicator
    oo_.([TYPE, 'TheoreticalMoments']).dsge.VarianceDecompositionME.density.(var).(exo) = density;
end