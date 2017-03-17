function indmcf = mcf_analysis(lpmat, ibeha, inobeha, options_mcf, DynareOptions)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu
%

% Copyright (C) 2014 European Commission
% Copyright (C) 2016 Dynare Team
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

pvalue_ks = options_mcf.pvalue_ks;
pvalue_corr = options_mcf.pvalue_corr;
alpha2 = options_mcf.alpha2;
param_names = options_mcf.param_names;

if DynareOptions.TeX
    if ~isfield(options_mcf,'param_names_tex')
        param_names_tex = options_mcf.param_names;
    else
        param_names_tex = options_mcf.param_names_tex;
    end
end
amcf_name = options_mcf.amcf_name;
amcf_title = options_mcf.amcf_title;
beha_title = options_mcf.beha_title;
nobeha_title = options_mcf.nobeha_title;
title = options_mcf.title;
fname_ = options_mcf.fname_;
xparam1=[];
if isfield(options_mcf,'xparam1'),
    xparam1=options_mcf.xparam1;
end    
OutputDirectoryName = options_mcf.OutputDirectoryName;

[proba, dproba] = stab_map_1(lpmat, ibeha, inobeha, [],0);
%         indindet=find(dproba>ksstat);
indmcf=find(proba<pvalue_ks);
[tmp,jtmp] = sort(proba(indmcf),2,'ascend');
indmcf = indmcf(jtmp);
if ~isempty(indmcf)
    skipline()
    headers=char('Parameter','d-stat','p-value');
    labels=char(param_names(indmcf,:));
    data_mat=[dproba(indmcf)' proba(indmcf)'];
    options_temp.noprint=0;
    dyntable(options_temp,['Smirnov statistics in driving ', title],headers,labels,data_mat,size(labels,2)+2,16,3);
    if DynareOptions.TeX
        labels_TeX=param_names_tex(indmcf,:);
        M_temp.dname=OutputDirectoryName ;
        M_temp.fname=fname_;
        dyn_latex_table(M_temp,options_temp,['Smirnov statistics in driving ', strrep(title,'_','\\_')],amcf_name,headers,labels_TeX,data_mat,size(labels,2)+2,16,6);
    end
end
    

if length(ibeha)>10 && length(inobeha)>10,
    indcorr1 = stab_map_2(lpmat(ibeha,:),alpha2, pvalue_corr, beha_title);
    indcorr2 = stab_map_2(lpmat(inobeha,:),alpha2, pvalue_corr, nobeha_title);
    indcorr = union(indcorr1(:), indcorr2(:));
    indcorr = indcorr(~ismember(indcorr(:),indmcf));
    indmcf = [indmcf(:); indcorr(:)];
end
if ~isempty(indmcf) && ~DynareOptions.nograph,
    skipline()
    xx=[];
    if ~ isempty(xparam1), xx=xparam1(indmcf); end
    scatter_mcf(lpmat(ibeha,indmcf),lpmat(inobeha,indmcf), param_names(indmcf,:), ...
        '.', [fname_,'_',amcf_name], OutputDirectoryName, amcf_title,xx, DynareOptions, ...
        beha_title, nobeha_title)
end
