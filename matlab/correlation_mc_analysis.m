function oo_ = correlation_mc_analysis(SampleSize,type,dname,fname,vartan,nvar,var1,var2,nar,mh_conf_sig,oo_,M_,options_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variables correlation function.

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

if isfield(oo_,[TYPE 'TheoreticalMoments'])
    eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments;'])
    if isfield(temporary_structure,'dsge')
        eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge;'])
        if isfield(temporary_structure,'correlation')
            eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge.correlation.Mean;'])
            if isfield(temporary_structure,deblank(var1))
                eval(['temporary_structure_1 = oo_.' TYPE 'TheoreticalMoments.dsge.correlation.Mean.' var1 ';']) 
                if isfield(temporary_structure_1,deblank(var2))
                    eval(['temporary_structure_2 = temporary_structure_1.' var2 ';'])
                    l1 = length(temporary_structure_2);
                    if l1<nar
                        % INITIALIZATION:
                        oo_ = initialize_output_structure(var1,var2,nar,type,oo_);
                        delete([PATH fname '_' TYPE 'Correlations*'])
                        [nvar,vartan,NumberOfFiles] = ...
                            dsge_simulated_theoretical_correlation(SampleSize,nar,M_,options_,oo_,type);
                    else
                        if ~isnan(temporary_structure_2(nar))
                            %Nothing to do.
                            return
                        end
                    end
                else
                    oo_ = initialize_output_structure(var1,var2,nar,TYPE,oo_);
                end
            else
                oo_ = initialize_output_structure(var1,var2,nar,TYPE,oo_);
            end
        else
            oo_ = initialize_output_structure(var1,var2,nar,TYPE,oo_);
        end
    else
        oo_ = initialize_output_structure(var1,var2,nar,TYPE,oo_);
    end
else
    oo_ = initialize_output_structure(var1,var2,nar,TYPE,oo_);
end
ListOfFiles = dir([ PATH  fname '_' TYPE 'Correlations*.mat']);
i1 = 1; tmp = zeros(SampleSize,1);
for file = 1:length(ListOfFiles)
    load([ PATH  ListOfFiles(file).name ]);
    i2 = i1 + rows(Correlation_array) - 1;
    tmp(i1:i2) = Correlation_array(:,indx1,indx2,nar);
    i1 = i2+1;
end
name = [ var1 '.' var2 ];
if ~isconst(tmp)
    [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
    if isfield(oo_,[ TYPE 'TheoreticalMoments'])
        eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments;'])
        if isfield(temporary_structure,'dsge')
            eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge;'])
            if isfield(temporary_structure,'correlation')
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Mean',nar,p_mean);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Median',nar,p_median);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Variance',nar,p_var);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'HPDinf',nar,hpd_interval(1));
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'HPDsup',nar,hpd_interval(2));
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'deciles',nar,p_deciles);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'density',nar,density);
            end
        end
    end
else
    if isfield(oo_,'PosteriorTheoreticalMoments')
        eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments;'])
        if isfield(temporary_structure,'dsge')
            eval(['temporary_structure = oo_.' TYPE 'TheoreticalMoments.dsge;'])
            if isfield(temporary_structure,'correlation')
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Mean',nar,NaN);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Median',nar,NaN);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Variance',nar,NaN);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'HPDinf',nar,NaN);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'HPDsup',nar,NaN);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'deciles',nar,NaN);
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'density',nar,NaN);
            end
        end
    end
end

function oo_ = initialize_output_structure(var1,var2,nar,type,oo_)
name = [ var1 '.' var2 ];
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.Mean.' name ' = NaN(' int2str(nar) ',1);']);
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.Median.' name ' = NaN(' int2str(nar) ',1);']);
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.Variance.' name ' = NaN(' int2str(nar) ',1);']);
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.HPDinf.' name ' = NaN(' int2str(nar) ',1);']);
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.HPDsup.' name ' = NaN(' int2str(nar) ',1);']);
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.deciles.' name ' = cell(' int2str(nar) ',1);']);
eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.density.' name ' = cell(' int2str(nar) ',1);']);
for i=1:nar
    eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.density.' name '(' int2str(i) ',1) = {NaN};']);
    eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.deciles.' name '(' int2str(i) ',1) = {NaN};']);
end

function oo_ = fill_output_structure(var1,var2,type,oo_,moment,lag,result)
name = [ var1 '.' var2 ];
switch moment
  case {'Mean','Median','Variance','HPDinf','HPDsup'} 
    eval(['oo_.' type  'TheoreticalMoments.dsge.correlation.' moment '.' name '(' int2str(lag) ',1) = result;']);
  case {'deciles','density'}
    eval(['oo_.' type 'TheoreticalMoments.dsge.correlation.' moment '.' name '(' int2str(lag) ',1) = {result};']);
  otherwise
    disp('fill_output_structure:: Unknown field!')
end