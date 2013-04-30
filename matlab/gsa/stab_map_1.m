function [proba, dproba] = stab_map_1(lpmat, ibehaviour, inonbehaviour, aname, iplot, ipar, dirname, pcrit, atitle)
%function [proba, dproba] = stab_map_1(lpmat, ibehaviour, inonbehaviour, aname, iplot, ipar, dirname, pcrit)
%
% lpmat =  Monte Carlo matrix
% ibehaviour = index of behavioural runs
% inonbehaviour = index of non-behavioural runs
% aname = label of the analysis
% iplot = 1 plot cumulative distributions (default)
% iplot = 0 no plots
% ipar = index array of parameters to plot
% dirname = (OPTIONAL) path of the directory where to save 
%            (default: current directory)
% pcrit = (OPTIONAL) critical value of the pvalue below which show the plots
%
% Plots: dotted lines for BEHAVIOURAL
%        solid lines for NON BEHAVIOURAL
% USES smirnov
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

global estim_params_ bayestopt_ M_ options_

if nargin<5,
  iplot=1;
end
fname_ = M_.fname;
if nargin<7,
  dirname='';
end
if nargin<9,
  atitle=aname;
end

nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

npar=size(lpmat,2);
ishock= npar>estim_params_.np;

if nargin<6,
  ipar=[];
end
if nargin<8 || isempty(pcrit),
  pcrit=1;
end

% Smirnov test for Blanchard; 
for j=1:npar,
  [H,P,KSSTAT] = smirnov(lpmat(ibehaviour,j),lpmat(inonbehaviour,j));
  proba(j)=P;
  dproba(j)=KSSTAT;
end
if isempty(ipar),
%     ipar=find(dproba>dcrit);
    ipar=find(proba<pcrit);
end
nparplot=length(ipar);
if iplot && ~options_.nograph
  lpmat=lpmat(:,ipar);
  ftit=bayestopt_.name(ipar+nshock*(1-ishock));
  
for i=1:ceil(nparplot/12),
  hh=dyn_figure(options_,'name',atitle);
  for j=1+12*(i-1):min(nparplot,12*i),
    subplot(3,4,j-12*(i-1))
    if ~isempty(ibehaviour),
      h=cumplot(lpmat(ibehaviour,j));
      set(h,'color',[0 0 1], 'linestyle',':','LineWidth',1.5)
    end
    hold on,
    if ~isempty(inonbehaviour),
      h=cumplot(lpmat(inonbehaviour,j));
      set(h,'color',[0 0 0],'LineWidth',1.5)
    end
%     title([ftit{j},'. D-stat ', num2str(dproba(ipar(j)),2)],'interpreter','none')
    title([ftit{j},'. p-value ', num2str(proba(ipar(j)),2)],'interpreter','none')
  end
  dyn_saveas(hh,[dirname,filesep,fname_,'_',aname,'_SA_',int2str(i)],options_);
end
end
