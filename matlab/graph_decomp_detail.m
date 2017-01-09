function []=graph_decomp_detail(z,shock_names,endo_names,i_var,initial_date,DynareModel,DynareOptions,opts_decomp)
%function []=graph_decomp_detail(z,shock_names,endo_names,i_var,initial_date,DynareModel,DynareOptions)
% Plots the results from the shock_decomposition command
% 
% Inputs
%   z               [n_var*(nshock+2)*nperiods]     shock decomposition array, see shock_decomposition.m for details
%   shock_names     [endo_nbr*string length]        shock names from M_.exo_names
%   endo_names      [exo_nbr*string length]         variable names from M_.endo_names
%   i_var           [n_var*1]                       vector indices of requested variables in M_.endo_names and z
%   initial_date    [dseries object]                first period of decomposition to plot
%   DynareModel     [structure]                     Dynare model structure
%   DynareOptions   [structure]                     Dynare options structure

% Copyright (C) 2010-2016 Dynare Team
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

% interactive = 0;
fig_mode='';
fig_mode1='';
% fig_names='';
% screen_shocks=0;
use_shock_groups = DynareOptions.use_shock_groups;
if use_shock_groups
    shock_groups = DynareModel.shock_groups.(use_shock_groups);
    shock_ind = fieldnames(shock_groups);
end

% number of components equals number of shocks + 1 (initial conditions)
comp_nbr = size(z,2)-1;

opts_decomp = DynareOptions.shock_decomp;

interactive = opts_decomp.interactive;
if ~isempty(opts_decomp.fig_mode)
    fig_mode = opts_decomp.fig_mode;
    fig_mode1 = ['_' fig_mode];
    fig_mode = [fig_mode '_'];
end
if DynareOptions.use_shock_groups,
    screen_shocks=0;
elseif comp_nbr>18
    screen_shocks = opts_decomp.screen_shocks;
end
fig_names = opts_decomp.fig_names;
%         fig_names = ['_' fig_names];
fig_names1 = [fig_names];
fig_names = [fig_names '_'];
if screen_shocks
    fig_names1 = [fig_names1 '_screen'];
    fig_names = [fig_names 'screen_'];
end


gend = size(z,3);
if isempty(initial_date)
    x = 0:gend;
    freq = 1;
else
    freq = initial_date.freq;
    initial_period = initial_date.time(1) + (initial_date.time(2)-1)/freq;
    x = initial_period-1/freq:(1/freq):initial_period+(gend-1)/freq;
end

ind_yrs = find(floor(x)==x); 
dind_tick = 1;
if floor(length(ind_yrs)/3);
    dind_tick = floor(length(ind_yrs)/3);
    xind_tick = x(ind_yrs(1)):dind_tick:x(ind_yrs(end))+(length(ind_yrs)-(dind_tick*3+1));
else
    xind_tick = x(ind_yrs(1)):dind_tick:x(ind_yrs(end))+(length(ind_yrs)-(dind_tick+1));
end
% xind_tick = floor(x(1))-floor(dind_tick/2):dind_tick:ceil(x(end))+ceil(dind_tick/2);
if abs(floor(x(1))-xind_tick(1))-abs(ceil(x(end))-xind_tick(end))>1
    xind_tick = xind_tick-1;
end
if abs(floor(x(1))-xind_tick(1))-abs(ceil(x(end))-xind_tick(end))<-1
    xind_tick = xind_tick+1;
end
% xind_tick = [x(ind_yrs(1))-floor(dind_tick/2):dind_tick:x(ind_yrs(end))+floor(dind_tick/2)]+1;
% xind_tick = x(ind_yrs(1))-1:dind_tick:x(ind_yrs(end))+1;
% xind_tick = x(ind_yrs(1))-1:dind_tick:x(ind_yrs(end))+dind_tick;

nvar = length(i_var);

%% write LaTeX-Header
if DynareOptions.TeX && any(strcmp('eps',cellstr(DynareOptions.graph_format)))
    fidTeX = fopen([DynareModel.fname '_shock_decomp' fig_mode1 fig_names1 '_detail.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by Dynare''s graph_decomp_detail.m.\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

ncol=3;
nrow=ceil(comp_nbr/ncol);
ntotrow = nrow;
nrow = min(ntotrow, 6);
nfigs = ceil(ntotrow/nrow);
labels = char(char(shock_names),'Initial values');
if ~(screen_shocks && comp_nbr>18),
    screen_shocks=0;
end
comp_nbr0=comp_nbr;
%%plot decomposition
for j=1:nvar
    z1 = squeeze(z(i_var(j),:,:));
    if screen_shocks,
        [junk, isort] = sort(mean(abs(z1(1:end-2,:)')), 'descend');
        labels = char(char(shock_names(isort(1:16),:)),'Others', 'Initial values');
        zres = sum(z1(isort(17:end),:),1);
        z1 = [z1(isort(1:16),:); zres; z1(comp_nbr0:end,:)];
        comp_nbr=18;
        nfigs=1;
    end
    xmin = x(1);
    xmin = min(xmin, xind_tick(1));
    xmax = x(end)+1/freq;
    xmax = max(xmax, xind_tick(end));
    ix = z1(1:comp_nbr,:) > 0;
    ymax = max(sum(z1(1:comp_nbr,:).*ix))*1.1;
    ix = z1(1:comp_nbr,:) < 0;
    ymin = min(sum(z1(1:comp_nbr,:).*ix))*1.1;
    if ymax-ymin < 1e-6
        continue
    end
    for jf = 1:nfigs
    fhandle = dyn_figure(DynareOptions,'Name',['Shock decomposition (detail): ' endo_names(i_var(j),:) fig_mode fig_names1],'position',[200 100 650 850], 'PaperPositionMode', 'auto','PaperOrientation','portrait','renderermode','auto');
    a0=zeros(1,4);
    a0(3)=inf;
    a0(4)=-inf;
    for ic=1+nrow*ncol*(jf-1):min(nrow*ncol*jf,comp_nbr),
        i = ic-nrow*ncol*(jf-1);
        zz = z1(ic,:);        
        zz(2,:)=z1(end,:)-zz;
        ipos=zz>0;
        ineg=zz<0;
        hax = subplot(nrow,ncol,i); set(gca,'box','on')
        hbar = bar(x(2:end),(zz.*ipos)','stacked');
        colormap([0.15 0.15 0.15;0.85 0.85 0.85]),
        set(hbar,'edgecolor','flat');
        hold on,
        hbar = bar(x(2:end),(zz.*ineg)','stacked');
        colormap([0.15 0.15 0.15;0.85 0.85 0.85]),
        set(hbar,'edgecolor','flat');
        title(deblank(labels(ic,:)),'Interpreter','none'),
        axis tight;
        a=axis;
        set(gca,'Xtick',xind_tick)
        set(gca,'xlim',[xmin xmax])
        a0(3)=min(a(3),a0(3));
        a0(4)=max(a(4),a0(4));
        set(gca,'ylim',a0(3:4))
        hold on, h1=plot(x(2:end),z1(end,:),'k-','LineWidth',2);
        if interactive & (~isoctave & use_shock_groups)
            mydata.use_shock_groups = DynareOptions.use_shock_groups;
            mydata.group_members = shock_groups.(shock_ind{ic}).shocks;
            set(gca,'UserData',mydata);
            if ~isempty(mydata.group_members{1})
                c = uicontextmenu;
                hax.UIContextMenu=c;
                browse_menu = uimenu(c,'Label','Browse group');
                expand_menu = uimenu(c,'Label','Expand group','Callback',['expand_group(''' mydata.use_shock_groups ''',''' deblank(endo_names(i_var(j),:)) ''',' int2str(ic) ')']);
                for jmember = mydata.group_members
                    uimenu('parent',browse_menu,'Label',char(jmember))
                end
            end
        end
    end
    for isub=1:i,
        subplot(nrow,ncol,isub),
        set(gca,'ylim',a0(3:4))
    end
    if nfigs>1,
        suffix = ['detail_' int2str(jf)];
    else
        suffix = ['detail'];
    end
    dyn_saveas(fhandle,[DynareModel.fname,'_shock_decomposition_',fig_mode,deblank(endo_names(i_var(j),:)),fig_names suffix],DynareOptions);
    if DynareOptions.TeX && any(strcmp('eps',cellstr(DynareOptions.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,['\\includegraphics[width=0.8\\textwidth]{%s_shock_decomposition_%s}\n'],DynareModel.fname,[fig_mode deblank(endo_names(i_var(j),:)) fig_names suffix]);
        fprintf(fidTeX,'\\label{Fig:shock_decomp:%s}\n',deblank(endo_names(i_var(j),:)));
        fprintf(fidTeX,'\\caption{Historical shock decomposition: $ %s $}\n',deblank(DynareModel.endo_names_tex(i_var(j),:)));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end    
    end
end

%% write LaTeX-Footer
if DynareOptions.TeX && any(strcmp('eps',cellstr(DynareOptions.graph_format)))
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
end
