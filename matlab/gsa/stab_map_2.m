function stab_map_2(x,alpha2, pvalue, fnam, dirname,xparam1,figtitle)
% function stab_map_2(x, alpha2, pvalue, fnam, dirname,xparam1)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.

% Copyright (C) 2012 Dynare Team
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

%global bayestopt_ estim_params_ dr_ options_ ys_ fname_
global bayestopt_ estim_params_ options_ oo_ M_

npar=size(x,2);
nsam=size(x,1);
ishock= npar>estim_params_.np;
if nargin<4,
  fnam='';
end
if nargin<5,
  dirname='';
end
if nargin<6,
  xparam1=[];
end
if nargin<7,
  figtitle=fnam;
end

ys_ = oo_.dr.ys;
dr_ = oo_.dr;
fname_ = M_.fname;
nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

c0=corrcoef(x);
c00=tril(c0,-1);
fig_nam_=[fname_,'_',fnam,'_corr_'];

ifig=0;
j2=0;
if ishock==0
  npar=estim_params_.np;
  if ~isempty(xparam1),
      xparam1=xparam1(nshock+1:end);
  end
else
  npar=estim_params_.np+nshock;
end
disp([' '])
disp(['Correlation analysis for ',fnam])

for j=1:npar,
    i2=find(abs(c00(:,j))>alpha2);
    if length(i2)>0,
        for jx=1:length(i2),
            tval  = abs(c00(i2(jx),j)*sqrt( (nsam-2)/(1-c00(i2(jx),j)^2) ));
            tcr = tcrit(nsam-2,pvalue);
            if tval>tcr,
                j2=j2+1;
                if ishock,
                    tmp_name = (['[',bayestopt_.name{j},',',bayestopt_.name{i2(jx)},']']);
                else
                    tmp_name = (['[',bayestopt_.name{j+nshock},',',bayestopt_.name{i2(jx)+nshock},']']);
                end
                fprintf(1,'%20s: corrcoef = %7.3f\n',tmp_name,c0(i2(jx),j));
                    
                if ~options_.nograph,
                    
                if mod(j2,12)==1,
                    ifig=ifig+1;
                    hh=dyn_figure(options_,'name',['Correlations in the ',figtitle,' sample ', num2str(ifig)]);
                end
                subplot(3,4,j2-(ifig-1)*12)
                %             bar(c0(i2,j)),
                %             set(gca,'xticklabel',bayestopt_.name(i2)),
                %             set(gca,'xtick',[1:length(i2)])
                %plot(stock_par(ixx(nfilt+1:end,i),j),stock_par(ixx(nfilt+1:end,i),i2(jx)),'.k')
                %hold on,
                plot(x(:,j),x(:,i2(jx)),'.')
                if ~isempty(xparam1)
                    hold on, plot(xparam1(j),xparam1(i2(jx)),'ro')
                end
                %             xlabel(deblank(estim_params_.param_names(j,:)),'interpreter','none'),
                %             ylabel(deblank(estim_params_.param_names(i2(jx),:)),'interpreter','none'),
                if ishock,
                    xlabel(bayestopt_.name{j},'interpreter','none'),
                    ylabel(bayestopt_.name{i2(jx)},'interpreter','none'),
                else
                    xlabel(bayestopt_.name{j+nshock},'interpreter','none'),
                    ylabel(bayestopt_.name{i2(jx)+nshock},'interpreter','none'),
                end
                title(['cc = ',num2str(c0(i2(jx),j))])
                if (mod(j2,12)==0) && j2>0,
                    dyn_saveas(hh,[dirname,filesep,fig_nam_,int2str(ifig)],options_);
                end
                end
            end
            
        end
    end
    if ~options_.nograph && (j==(npar)) && j2>0 && (mod(j2,12)~=0),
        dyn_saveas(hh,[dirname,filesep,fig_nam_,int2str(ifig)],options_);
    end
    
end
if j2==0,
    disp(['No correlation term with pvalue <', num2str(pvalue),' found for ',fnam])
end
%close all
