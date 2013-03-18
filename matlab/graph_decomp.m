function []=graph_decomp(z,shock_names,endo_names,i_var,initial_date,DynareModel,DynareOptions)
%function []=graph_decomp(z,varlist,initial_period,freq)

% Copyright (C) 2010-2011 Dynare Team
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

% number of components equals number of shocks + 1 (initial conditions)
comp_nbr = size(z,2)-1;

gend = size(z,3);
freq = initial_date.freq;
initial_period = initial_date.period + initial_date.subperiod/freq;
x = initial_period-1/freq:(1/freq):initial_period+(gend-1)/freq;

nvar = length(i_var);

for j=1:nvar
    z1 = squeeze(z(i_var(j),:,:));
    xmin = x(1);
    xmax = x(end);
    ix = z1 > 0;
    ymax = max(sum(z1.*ix));
    ix = z1 < 0;
    ymin = min(sum(z1.*ix));
    if ymax-ymin < 1e-6
        continue
    end
    fhandle = dyn_figure(DynareOptions,'Name',['Shock decomposition: ',endo_names(i_var(j),:)]);
    ax=axes('Position',[0.1 0.1 0.6 0.8]);
    axis(ax,[xmin xmax ymin ymax]);
    plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
    hold on;
    for i=1:gend
        i_1 = i-1;
        yp = 0;
        ym = 0;
        for k = 1:comp_nbr
            zz = z1(k,i);
            if zz > 0
                fill([x(i) x(i) x(i+1) x(i+1)],[yp yp+zz yp+zz yp],k);
                yp = yp+zz;
            else
                fill([x(i) x(i) x(i+1) x(i+1)],[ym ym+zz ym+zz ym],k);
                ym = ym+zz;
            end
            hold on;
        end
    end
    plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
    hold off;

    axes('Position',[0.75 0.1 0.2 0.8]);
    axis([0 1 0 1]);
    axis off;
    hold on;
    y1 = 0;
    height = 1/comp_nbr;
    labels = char(shock_names,'Initial values');

    for i=1:comp_nbr
        fill([0 0 0.2 0.2],[y1 y1+0.7*height y1+0.7*height y1],i);
        hold on
        text(0.3,y1+0.3*height,labels(i,:),'Interpreter','none');
        hold on
        y1 = y1 + height;
    end

    dyn_saveas(fhandle,[DynareModel.fname,'_shock_decomposition_',deblank(endo_names(i_var(j),:))],DynareOptions);
    hold off
end