function pick
%
% Copyright (C) 2001-2017 European Commission
% 
% This file is part of GLUEWIN
% GLUEWIN is a MATLAB code designed for analysing the output
% of Monte Carlo runs when empirical observations of the model output are available
% and implements the GSA-GLUE methodology by Ratto et al. [1], based on a combination
% of GLUE (Generalised Likelihood Uncertainty Estimation) by K. Beven [2] and GSA
% Global Sensitivity Analysis) [3].']
% The program has been developed by M. Ratto, European Commission, Joint Research Centre,
% Institute for the Protection and Security of The Citizen, Technological and Economic Risk Management,
% Applied Statistics, as a deliverable of the IMPACT project 
% (EC Fifth Framework Programme, SCA Project, IST-1999-11313, DG-INFSO).
%
% The graphical layout of the code is inspired by the freeware GLUE package by K. Beven,
% vailable at the Lancaster University web site on the page [4]:
% http://www.es.lancs.ac.uk/hfdg/glue.html
% to which the GLUEWIN code introduces several extensions and additional options.
% Thanks are due to R. Girardi, A. Rossi, A. Saltelli, S. Tarantola and U. Callies for numerous
% comments and suggestions.
% For more information, please contact marco.ratto@ec.europa.eu
% 
% Disclaimer: This software has been developed at the Joint Research Centre of European Commission
% by officers in the course of their official duties. This software is not subject to copyright
% protection and is in the public domain. It is an experimental system. The Joint Research Centre
% of European Commission assumes no responsibility whatsoever for its use by other parties
% and makes no guarantees, expressed or implied, about its quality, reliability, or any other
% characteristic. We would appreciate acknowledgement if the software is used.
%
% [1] Ratto, M., Tarantola, S., A. Saltelli, Sensitivity analysis in model calibration: GSA-GLUE approach.
%                'Computer Physics Communications, 136, 2001, 212-224
% [2] Beven K.J., Binley A., The Future of Distributed Models: Model Calibration and Uncertainty
%                'Prediction, Hydrological Processes, 6, 279-298, 1992
% [3] Saltelli, A., K. Chan, M. Scott, Editors, (2000), Sensitivity analysis, John Wiley & Sons
%                'publishers, Probability and Statistics series.
% [4] Beven K., GLUE for Windows User manual, 1998.      



pmenu=findobj(gcf,'type','uicontextmenu','Tag','Run viewer');
button1=findobj(gcf,'type','uimenu','Tag','save params');
button2=findobj(gcf,'type','uimenu','Tag','eval params');
%button=get(pmenu,'children');
gg=gco;
ax0=gca;
set(gg,'buttondownfcn',[]);
c=get(gca,'currentpoint');
x=c(1,1);
y=c(1,2);
X=get(gco,'xdata');
Y=get(gco,'ydata');
dx=get(gca,'xlim');
dy=get(gca,'ylim');
pos=get(gca,'position');
scalex=dx(2)-dx(1);
scaley=dy(2)-dy(1);
if length(X)>1,
    K = dsearchn([(Y./scaley)' (X./scalex)'],[y/scaley x/scalex]);
else
    az=get(gca,'children');
    T =get(az(end),'ydata');
    [dum K]=max(T);
end

KK=K;

set(button1,'Label',['Save ',num2str(K)],'Callback',['scatter_callback(',num2str(KK),',''save'')']);
set(button2,'Label',['Eval ',num2str(K)],'Callback',['scatter_callback(',num2str(KK),',''eval'')']);
hh=findobj(gcf,'type','axes','Tag','scatter');
for k=1:length(hh),
    axes(hh(k));
    dum=get(gca,'children');
    dumx=get(dum(end),'xdata');
    dumy=get(dum(end),'ydata');
    xmid=min(dumx) + 0.5*(max(dumx)-min(dumx));
    hold on
    plot(dumx(KK),dumy(KK),'or');     
    if dumx(KK) < xmid,
        text(dumx(KK),dumy(KK),['  ',num2str(K)], ...
            'FontWeight','Bold',...
            'Color','r');
    else
        text(dumx(KK),dumy(KK),[num2str(K),'  '], ...
            'HorizontalAlignment','right', ...
            'FontWeight','Bold',...
            'Color','r');
    end
    hold off
end