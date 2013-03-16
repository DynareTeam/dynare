function plot_identification(params,idemoments,idehess,idemodel, idelre, advanced, tittxt, name, IdentifDirectoryName)
% function plot_identification(params,idemoments,idehess,idemodel, idelre, advanced, tittxt, name, IdentifDirectoryName)
%
% INPUTS
%    o params             [array] parameter values for identification checks
%    o idemoments         [structure] identification results for the moments
%    o idehess            [structure] identification results for the Hessian
%    o idemodel           [structure] identification results for the reduced form solution
%    o idelre             [structure] identification results for the LRE model
%    o advanced           [integer] flag for advanced identification checks
%    o tittxt             [char] name of the results to plot 
%    o name               [char] list of names
%    o IdentifDirectoryName   [char] directory name
%    
% OUTPUTS
%    None
%    
% SPECIAL REQUIREMENTS
%    None

% Copyright (C) 2008-2012 Dynare Team
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

global M_ options_


[SampleSize, nparam]=size(params);
siJnorm = idemoments.siJnorm;
siHnorm = idemodel.siHnorm;
siLREnorm = idelre.siLREnorm;

% if prior_exist,
%     tittxt = 'Prior mean - ';
% else
%     tittxt = '';
% end
tittxt1=regexprep(tittxt, ' ', '_');
tittxt1=strrep(tittxt1, '.', '');
if SampleSize == 1,
    siJ = idemoments.siJ;
    hh = dyn_figure(options_,'Name',[tittxt, ' - Identification using info from observables']);
    subplot(211)
    mmm = (idehess.ide_strength_J);
    [ss, is] = sort(mmm);
    bar(log([idehess.ide_strength_J(:,is)' idehess.ide_strength_J_prior(:,is)']))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    legend('relative to param value','relative to prior std','Location','Best')
    if  idehess.flag_score,
        title('Identification strength with asymptotic Information matrix (log-scale)')
    else
        title('Identification strength with moments Information matrix (log-scale)')
    end
    
    subplot(212)
    bar(log([idehess.deltaM(is) idehess.deltaM_prior(is)]))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    legend('relative to param value','relative to prior std','Location','Best')
    if  idehess.flag_score,
        title('Sensitivity component with asymptotic Information matrix (log-scale)')
    else
        title('Sensitivity component with moments Information matrix (log-scale)')
    end
    dyn_saveas(hh,[IdentifDirectoryName '/' M_.fname '_ident_strength_' tittxt1],options_);
    
    if advanced,
        if ~options_.nodisplay,
            disp(' ')
            disp('Press ENTER to plot advanced diagnostics'), pause(5),
        end
        hh = dyn_figure(options_,'Name',[tittxt, ' - Sensitivity plot']);
        subplot(211)
        mmm = (siJnorm)'./max(siJnorm);
        mmm1 = (siHnorm)'./max(siHnorm);
        mmm=[mmm mmm1];
        mmm1 = (siLREnorm)'./max(siLREnorm);
        offset=length(siHnorm)-length(siLREnorm);
        mmm1 = [NaN(offset,1); mmm1];
        mmm=[mmm mmm1];
        
        bar(log(mmm(is,:).*100))
        set(gca,'xlim',[0 nparam+1])
        set(gca,'xticklabel','')
        dy = get(gca,'ylim');
        for ip=1:nparam,
            text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        legend('Moments','Model','LRE model','Location','Best')
        title('Sensitivity bars using derivatives (log-scale)')
        dyn_saveas(hh,[IdentifDirectoryName '/' M_.fname '_sensitivity_' tittxt1 ],options_);
        
        % identificaton patterns
        for  j=1:size(idemoments.cosnJ,2),
            pax=NaN(nparam,nparam);
%             fprintf('\n')
%             disp(['Collinearity patterns with ', int2str(j) ,' parameter(s)'])
%             fprintf('%-15s [%-*s] %10s\n','Parameter',(15+1)*j,' Expl. params ','cosn')
            for i=1:nparam,
                namx='';
                for in=1:j,
                    dumpindx = idemoments.pars{i,j}(in);
                    if isnan(dumpindx),
                        namx=[namx ' ' sprintf('%-15s','--')];
                    else
                        namx=[namx ' ' sprintf('%-15s',name{dumpindx})];
                        pax(i,dumpindx)=idemoments.cosnJ(i,j);
                    end
                end
%                 fprintf('%-15s [%s] %10.3f\n',name{i},namx,idemoments.cosnJ(i,j))
            end
            hh = dyn_figure(options_,'Name',[tittxt,' - Collinearity patterns with ', int2str(j) ,' parameter(s)']);
            imagesc(pax,[0 1]);
            set(gca,'xticklabel','')
            set(gca,'yticklabel','')
            for ip=1:nparam,
                text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
                text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
            end
            colorbar;
            ax=colormap;
            ax(1,:)=[0.9 0.9 0.9];
            colormap(ax);
            if nparam>10,
                set(gca,'xtick',(5:5:nparam))
                set(gca,'ytick',(5:5:nparam))
            end
            set(gca,'xgrid','on')
            set(gca,'ygrid','on')
            xlabel([tittxt,' - Collinearity patterns with ', int2str(j) ,' parameter(s)'],'interpreter','none')
            dyn_saveas(hh,[ IdentifDirectoryName '/' M_.fname '_ident_collinearity_' tittxt1 '_' int2str(j) ],options_);
        end
        disp('')
        if idehess.flag_score,
            [U,S,V]=svd(idehess.AHess,0);
            S=diag(S);
            if nparam<5,
                f1 = dyn_figure(options_,'Name',[tittxt,' - Identification patterns (Information matrix)']);
            else
                f1 = dyn_figure(options_,'Name',[tittxt,' - Identification patterns (Information matrix): SMALLEST SV']);
                f2 = dyn_figure(options_,'Name',[tittxt,' - Identification patterns (Information matrix): HIGHEST SV']);
            end
        else
            S = idemoments.S;
            V = idemoments.V;
            if nparam<5,
                f1 = dyn_figure(options_,'Name',[tittxt,' - Identification patterns (moments)']);
            else
                f1 = dyn_figure(options_,'Name',[tittxt,' - Identification patterns (moments): SMALLEST SV']);
                f2 = dyn_figure(options_,'Name',[tittxt,' - Identification patterns (moments): HIGHEST SV']);
            end
        end
        for j=1:min(nparam,8),
            if j<5,
                set(0,'CurrentFigure',f1),
                jj=j;
            else
                set(0,'CurrentFigure',f2),
                jj=j-4;
            end
            subplot(4,1,jj),
            if j<5
                bar(abs(V(:,end-j+1))),
                Stit = S(end-j+1);
            else
                bar(abs(V(:,jj))),
                Stit = S(jj);
            end
            set(gca,'xticklabel','')
            if j==4 || j==nparam || j==8,
                for ip=1:nparam,
                    text(ip,-0.02,name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                end
            end
            title(['Singular value ',num2str(Stit)])
        end
        dyn_saveas(f1,[  IdentifDirectoryName '/' M_.fname '_ident_pattern_' tittxt1 '_1' ],options_);
        if nparam>4,
            dyn_saveas(f2,[  IdentifDirectoryName '/' M_.fname '_ident_pattern_' tittxt1 '_2' ],options_);
        end
    end
    
else
    hh = dyn_figure(options_,'Name',['MC sensitivities']);
    subplot(211)
    mmm = (idehess.ide_strength_J);
    [ss, is] = sort(mmm);
    mmm = mean(siJnorm)';
    mmm = mmm./max(mmm);
    if advanced,
        mmm1 = mean(siHnorm)';
        mmm=[mmm mmm1./max(mmm1)];
        mmm1 = mean(siLREnorm)';
        offset=size(siHnorm,2)-size(siLREnorm,2);
        mmm1 = [NaN(offset,1); mmm1./max(mmm1)];
        mmm=[mmm mmm1];
    end        
        
    bar(mmm(is,:))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    if advanced,
        legend('Moments','Model','LRE model','Location','Best')
    end
    title('MC mean of sensitivity measures')
    dyn_saveas(hh,[ IdentifDirectoryName '/' M_.fname '_MC_sensitivity' ],options_);
    if advanced,
        if ~options_.nodisplay,
            disp(' ')
            disp('Press ENTER to display advanced diagnostics'), pause(5),
        end
%         options_.nograph=1;
        hh = dyn_figure(options_,'Name','MC Condition Number');
        subplot(221)
        hist(log10(idemodel.cond))
        title('log10 of Condition number in the model')
        subplot(222)
        hist(log10(idemoments.cond))
        title('log10 of Condition number in the moments')
        subplot(223)
        hist(log10(idelre.cond))
        title('log10 of Condition number in the LRE model')
        dyn_saveas(hh,[IdentifDirectoryName '/' M_.fname '_ident_COND' ],options_);
        ncut=floor(SampleSize/10*9);
        [dum,is]=sort(idelre.cond);
        [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), 'MC_HighestCondNumberLRE', 1, [], IdentifDirectoryName, 0.1);
        [dum,is]=sort(idemodel.cond);
        [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), 'MC_HighestCondNumberModel', 1, [], IdentifDirectoryName, 0.1);
        [dum,is]=sort(idemoments.cond);
        [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), 'MC_HighestCondNumberMoments', 1, [], IdentifDirectoryName, 0.1);
%         [proba, dproba] = stab_map_1(idemoments.Mco', is(1:ncut), is(ncut+1:end), 'HighestCondNumberMoments_vs_Mco', 1, [], IdentifDirectoryName);
%         for j=1:nparam,
% %             ibeh=find(idemoments.Mco(j,:)<0.9);
% %             inonbeh=find(idemoments.Mco(j,:)>=0.9);
% %             if ~isempty(ibeh) && ~isempty(inonbeh)
% %                 [proba, dproba] = stab_map_1(params, ibeh, inonbeh, ['HighestMultiCollinearity_',name{j}], 1, [], IdentifDirectoryName);
% %             end
%             [~,is]=sort(idemoments.Mco(:,j));
%             [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), ['MC_HighestMultiCollinearity_',name{j}], 1, [], IdentifDirectoryName, 0.15);
%         end

        if nparam<5,
            f1 = dyn_figure(options_,'Name',[tittxt,' - MC Identification patterns (moments): HIGHEST SV']);
        else
            f1 = dyn_figure(options_,'Name',[tittxt,' - MC Identification patterns (moments): SMALLEST SV']);
            f2 = dyn_figure(options_,'Name',[tittxt,' - MC Identification patterns (moments): HIGHEST SV']);
        end
        nplots=min(nparam,8);
        if nplots>4,
            nsubplo=ceil(nplots/2);
        else
            nsubplo=nplots;
        end
        for j=1:nplots,
            if (nparam>4 && j<=ceil(nplots/2)) || nparam<5,
                set(0,'CurrentFigure',f1),
                jj=j;
                VVV=squeeze(abs(idemoments.V(:,:,end-j+1)));
                SSS = idemoments.S(:,end-j+1);
            else
                set(0,'CurrentFigure',f2),
                jj=j-ceil(nplots/2);
                VVV=squeeze(abs(idemoments.V(:,:,jj)));
                SSS = idemoments.S(:,jj);
            end
            subplot(nsubplo,1,jj),
            for i=1:nparam,
                [post_mean, post_median(:,i), post_var, hpd_interval(i,:), post_deciles] = posterior_moments(VVV(:,i),0,0.9);
            end
            bar(post_median)
            hold on, plot(hpd_interval,'--*r'),
            Stit=mean(SSS);

            set(gca,'xticklabel','')
            if j==4 || j==nparam || j==8,
                for ip=1:nparam,
                    text(ip,-0.02,name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                end
            end
            title(['MEAN Singular value ',num2str(Stit)])
        end
        dyn_saveas(f1,[IdentifDirectoryName '/' M_.fname '_MC_ident_pattern_1' ],options_);
        if nparam>4,
            dyn_saveas(f2,[  IdentifDirectoryName '/' M_.fname '_MC_ident_pattern_2' ],options_);
        end
    end
end

% disp_identification(params, idemodel, idemoments, name)