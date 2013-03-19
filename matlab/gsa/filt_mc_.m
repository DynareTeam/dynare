function [rmse_MC, ixx] = filt_mc_(OutDir,options_gsa_)
% function [rmse_MC, ixx] = filt_mc_(OutDir)
% inputs (from opt_gsa structure)
% vvarvecm = options_gsa_.var_rmse;
% loadSA   = options_gsa_.load_rmse;
% pfilt    = options_gsa_.pfilt_rmse;
% alpha    = options_gsa_.alpha_rmse;
% alpha2   = options_gsa_.alpha2_rmse;
% istart   = options_gsa_.istart_rmse;
% alphaPC  = 0.5;
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it
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

global bayestopt_ estim_params_ M_ options_ oo_

% options_gsa_=options_.opt_gsa;
vvarvecm = options_gsa_.var_rmse;
loadSA   = options_gsa_.load_rmse;
pfilt    = options_gsa_.pfilt_rmse;
alpha    = options_gsa_.alpha_rmse;
% alpha2   = options_gsa_.alpha2_rmse;
alpha2 = 0;
pvalue   = options_gsa_.alpha2_rmse;
istart   = options_gsa_.istart_rmse;
alphaPC  = 0.5;

fname_ = M_.fname;
lgy_ = M_.endo_names;
dr_ = oo_.dr;

disp(' ')
disp(' ')
disp('Starting sensitivity analysis')
disp('for the fit of EACH observed series ...')
disp(' ')
disp('Deleting old SA figures...')
a=dir([OutDir,filesep,'*.*']);
tmp1='0';
if options_.opt_gsa.ppost,
    tmp=['_rmse_post'];
else
    if options_.opt_gsa.pprior
        tmp=['_rmse_prior'];
    else
        tmp=['_rmse_mc'];
    end
    if options_gsa_.lik_only,
        tmp1 = [tmp,'_post_SA'];
        tmp = [tmp,'_lik_SA'];
    end
end
for j=1:length(a),
    if strmatch([fname_,tmp],a(j).name),
        disp(a(j).name)
        delete([OutDir,'/',a(j).name])
    end,
    if strmatch([fname_,tmp1],a(j).name),
        disp(a(j).name)
        delete([OutDir,'/',a(j).name])
    end,
end
disp('done !')


nshock=estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn;
npar=estim_params_.np;
if ~isempty(options_.mode_file),
    load(options_.mode_file,'xparam1'),
end
if options_.opt_gsa.ppost,
    c=load([fname_,'_mean.mat'],'xparam1');
    xparam1_mean=c.xparam1;
    clear c
elseif ~isempty(options_.mode_file) && exist([fname_,'_mean.mat'])==2,
    c=load([fname_,'_mean.mat'],'xparam1');
    xparam1_mean=c.xparam1;
    clear c
end

if options_.opt_gsa.ppost,
    fnamtmp=[fname_,'_post'];
    DirectoryName = CheckPath('metropolis',M_.dname);
else
    if options_.opt_gsa.pprior
        fnamtmp=[fname_,'_prior'];
        DirectoryName = CheckPath(['gsa' filesep 'prior'],M_.dname);
    else
        fnamtmp=[fname_,'_mc'];
        DirectoryName = CheckPath(['gsa' filesep 'mc'],M_.dname);
    end
end
if ~loadSA,
    if exist('xparam1','var')
        M_ = set_all_parameters(xparam1,estim_params_,M_);
        steady_(M_,options_,oo_);
        ys_mode=oo_.steady_state;
    end
    if exist('xparam1_mean','var')
        M_ = set_all_parameters(xparam1_mean,estim_params_,M_);
        steady_(M_,options_,oo_);
        ys_mean=oo_.steady_state;
    end
    %   eval(options_.datafile)
    obs = dat_fil_(options_.datafile);
    %stock_gend=data_info.gend;
    %stock_data = data_info.data;
    load([DirectoryName '/' M_.fname '_data.mat']);    
    filfilt = dir([DirectoryName filesep M_.fname '_filter_step_ahead*.mat']);
    filparam = dir([DirectoryName filesep M_.fname '_param*.mat']);
    x=[];
    logpo2=[];
    sto_ys=[];
    for j=1:length(filparam),
        %load([DirectoryName '/' M_.fname '_param',int2str(j),'.mat']);
        if isempty(strmatch([M_.fname '_param_irf'],filparam(j).name))
            load([DirectoryName filesep filparam(j).name]);
            x=[x; stock];
            logpo2=[logpo2; stock_logpo];
            sto_ys=[sto_ys; stock_ys];
            clear stock stock_logpo stock_ys;
        end
    end
    nruns=size(x,1);
    nfilt=floor(pfilt*nruns);
    if options_.opt_gsa.ppost || (options_.opt_gsa.ppost==0 && options_.opt_gsa.lik_only==0)
        disp(' ')
        disp('Computing RMSE''s...')
        fobs = options_.first_obs;
        nobs=options_.nobs;
        for i=1:size(vvarvecm,1),
            vj=deblank(vvarvecm(i,:));
            eval(['vobs =obs.',vj,'(fobs:fobs-1+nobs);'])
            if options_.prefilter == 1
                %eval([vj,'=',vj,'-bayestopt_.mean_varobs(i);'])
                %eval([vj,'=',vj,'-mean(',vj,',1);'])
                vobs = vobs-mean(vobs,1);
            end
            
            jxj = strmatch(vj,lgy_(dr_.order_var,:),'exact');
            js = strmatch(vj,lgy_,'exact');
            if exist('xparam1','var')
                %       if isfield(oo_,'FilteredVariables')
                %       eval(['rmse_mode(i) = sqrt(mean((vobs(istart:end)-oo_.steady_state(js)-oo_.FilteredVariables.',vj,'(istart:end-1)).^2));'])
                %       else
                [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(xparam1,stock_gend,stock_data,{},0);
                y0 = squeeze(aK(1,jxj,:)) + ...
                    kron(ys_mode(js,:),ones(size(aK,3),1));
                %       y0 = ahat(jxj,:)' + ...
                %         kron(ys_mode(js,:),ones(size(ahat,2),1));
                rmse_mode(i) = sqrt(mean((vobs(istart:end)-y0(istart:end-1)).^2));
                %       end
            end
            y0=zeros(nobs+1,nruns);
            nbb=0;
            for j=1:length(filfilt),
                load([DirectoryName filesep M_.fname '_filter_step_ahead',num2str(j),'.mat']);
                nb = size(stock,4);
                %         y0(:,nbb+1:nbb+nb)=squeeze(stock(1,js,:,:)) + ...
                %           kron(sto_ys(nbb+1:nbb+nb,js)',ones(size(stock,3),1));
                y0(:,nbb+1:nbb+nb)=squeeze(stock(1,js,1:nobs+1,:)) + ...
                    kron(sto_ys(nbb+1:nbb+nb,js)',ones(nobs+1,1));
                %y0(:,:,size(y0,3):size(y0,3)+size(stock,3))=stock;
                nbb=nbb+nb;
                clear stock;
            end
            y0M=mean(y0,2);
            for j=1:nruns,
                rmse_MC(j,i) = sqrt(mean((vobs(istart:end)-y0(istart:end-1,j)).^2));
            end
            if exist('xparam1_mean','var')
                %eval(['rmse_pmean(i) = sqrt(mean((',vj,'(fobs-1+istart:fobs-1+nobs)-y0M(istart:end-1)).^2));'])
                [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(xparam1_mean,stock_gend,stock_data,{},0);
                y0 = squeeze(aK(1,jxj,:)) + ...
                    kron(ys_mean(js,:),ones(size(aK,3),1));
                %       y0 = ahat(jxj,:)' + ...
                %         kron(ys_mean(js,:),ones(size(ahat,2),1));
                rmse_pmean(i) = sqrt(mean((vobs(istart:end)-y0(istart:end-1)).^2));
            end
        end
        clear stock_filter;
    end
    for j=1:nruns,
        lnprior(j,1) = priordens(x(j,:)',bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    end
    likelihood=logpo2(:)-lnprior(:);
    disp('... done!')
    
    if options_.opt_gsa.ppost
        save([OutDir,filesep,fnamtmp,'.mat'], 'x', 'logpo2', 'likelihood', 'rmse_MC', 'rmse_mode','rmse_pmean')
    else
        if options_.opt_gsa.lik_only
            save([OutDir,filesep,fnamtmp, '.mat'], 'x', 'logpo2','likelihood', '-append')
        else
            save([OutDir,filesep,fnamtmp, '.mat'], 'x', 'logpo2','likelihood', 'rmse_MC','-append')
            if exist('xparam1_mean','var')
                save([OutDir,filesep,fnamtmp, '.mat'], 'rmse_pmean','-append')
            end
            if exist('xparam1','var')
                save([OutDir,filesep,fnamtmp,'.mat'], 'rmse_mode','-append')
            end
        end
    end
else
    if options_.opt_gsa.lik_only && options_.opt_gsa.ppost==0
        load([OutDir,filesep,fnamtmp, '.mat'],'x','logpo2','likelihood');
    else
        load([OutDir,filesep,fnamtmp, '.mat'],'x','logpo2','likelihood','rmse_MC','rmse_mode','rmse_pmean');
    end
    lnprior=logpo2(:)-likelihood(:);
    nruns=size(x,1);
    nfilt=floor(pfilt*nruns);
end
% smirnov tests
nfilt0=nfilt*ones(size(vvarvecm,1),1);
logpo2=logpo2(:);
if ~options_.opt_gsa.ppost
    [dum, ipost]=sort(-logpo2);
    [dum, ilik]=sort(-likelihood);
end
if ~options_.opt_gsa.ppost && options_.opt_gsa.lik_only
    if options_.opt_gsa.pprior
        anam='rmse_prior_post';
    else
        anam='rmse_mc_post';
    end
    stab_map_1(x, ipost(1:nfilt), ipost(nfilt+1:end), anam, 1,[],OutDir);
    stab_map_2(x(ipost(1:nfilt),:),alpha2,pvalue,anam, OutDir);
    if options_.opt_gsa.pprior
        anam='rmse_prior_lik';
    else
        anam='rmse_mc_lik';
    end
    stab_map_1(x, ilik(1:nfilt), ilik(nfilt+1:end), anam, 1,[],OutDir);
    stab_map_2(x(ilik(1:nfilt),:),alpha2,pvalue,anam, OutDir);
else
    for i=1:size(vvarvecm,1),
        [dum, ixx(:,i)]=sort(rmse_MC(:,i));
        if options_.opt_gsa.ppost,
            %nfilt0(i)=length(find(rmse_MC(:,i)<rmse_pmean(i)));
            rmse_txt=rmse_pmean;
        else
            if options_.opt_gsa.pprior || ~exist('rmse_pmean'),
                if exist('rmse_mode'),
                    rmse_txt=rmse_mode;
                else
                    rmse_txt=NaN(1,size(rmse_MC,2));
                end
            else
                %nfilt0(i)=length(find(rmse_MC(:,i)<rmse_pmean(i)));
                rmse_txt=rmse_pmean;
            end
        end
        for j=1:npar+nshock,
            [H,P,KSSTAT] = smirnov(x(ixx(nfilt0(i)+1:end,i),j),x(ixx(1:nfilt0(i),i),j), alpha);
            [H1,P1,KSSTAT1] = smirnov(x(ixx(nfilt0(i)+1:end,i),j),x(ixx(1:nfilt0(i),i),j),alpha,1);
            [H2,P2,KSSTAT2] = smirnov(x(ixx(nfilt0(i)+1:end,i),j),x(ixx(1:nfilt0(i),i),j),alpha,-1);
            if H1 & H2==0,
                SS(j,i)=1;
            elseif H1==0,
                SS(j,i)=-1;
            else
                SS(j,i)=0;
            end
            PP(j,i)=P;
        end
    end
    if ~options_.nograph,
    ifig=0;
    for i=1:size(vvarvecm,1),
        if mod(i,9)==1,
            ifig=ifig+1;
            hh=dyn_figure(options_,'name',['Prior ',int2str(ifig)]);
        end
        subplot(3,3,i-9*(ifig-1))
        h=cumplot(lnprior(ixx(1:nfilt0(i),i)));
        set(h,'color','red')
        hold on, cumplot(lnprior)
        h=cumplot(lnprior(ixx(nfilt0(i)+1:end,i)));
        set(h,'color','green')
        title(vvarvecm(i,:),'interpreter','none')
        if mod(i,9)==0 || i==size(vvarvecm,1)
            if options_.opt_gsa.ppost
                dyn_saveas(hh,[OutDir '/' fname_ '_rmse_post_lnprior',int2str(ifig)],options_);
            else
                if options_.opt_gsa.pprior
                    dyn_saveas(hh,[OutDir '/' fname_ '_rmse_prior_lnprior',int2str(ifig) ],options_);
                else
                    dyn_saveas(hh,[OutDir '/' fname_ '_rmse_mc_lnprior',int2str(ifig) ],options_);
                end
            end
            if ~options_.nodisplay
                close(hh);
            end
        end
    end
    ifig=0;
    for i=1:size(vvarvecm,1),
        if mod(i,9)==1,
            ifig=ifig+1;
            hh = dyn_figure(options_,'Name',['Likelihood ',int2str(ifig)]);
        end
        subplot(3,3,i-9*(ifig-1))
        h=cumplot(likelihood(ixx(1:nfilt0(i),i)));
        set(h,'color','red')
        hold on, h=cumplot(likelihood);
        h=cumplot(likelihood(ixx(nfilt0(i)+1:end,i)));
        set(h,'color','green')
        title(vvarvecm(i,:),'interpreter','none')
        if options_.opt_gsa.ppost==0,
            set(gca,'xlim',[min( likelihood(ixx(1:nfilt0(i),i)) ) max( likelihood(ixx(1:nfilt0(i),i)) )])
        end
        if mod(i,9)==0 || i==size(vvarvecm,1)
            if options_.opt_gsa.ppost
                dyn_saveas(hh,[OutDir '/' fname_ '_rmse_post_lnlik',int2str(ifig) ],options_);
            else
                if options_.opt_gsa.pprior
                    dyn_saveas(hh,[OutDir '/' fname_ '_rmse_prior_lnlik',int2str(ifig)],options_);
                else
                    dyn_saveas(hh,[OutDir '/' fname_ '_rmse_mc_lnlik',int2str(ifig) ],options_);
                end
            end
            if ~options_.nodisplay
                close(hh);
            end
        end
    end
    ifig=0;
    for i=1:size(vvarvecm,1),
        if mod(i,9)==1,
            ifig=ifig+1;
            hh = dyn_figure(options_,'Name',['Posterior ',int2str(ifig)]);
        end
        subplot(3,3,i-9*(ifig-1))
        h=cumplot(logpo2(ixx(1:nfilt0(i),i)));
        set(h,'color','red')
        hold on, h=cumplot(logpo2);
        h=cumplot(logpo2(ixx(nfilt0(i)+1:end,i)));
        set(h,'color','green')
        title(vvarvecm(i,:),'interpreter','none')
        if options_.opt_gsa.ppost==0,
            set(gca,'xlim',[min( logpo2(ixx(1:nfilt0(i),i)) ) max( logpo2(ixx(1:nfilt0(i),i)) )])
        end
        if mod(i,9)==0 || i==size(vvarvecm,1)
            if options_.opt_gsa.ppost
                dyn_saveas(hh,[OutDir '/' fname_ '_rmse_post_lnpost',int2str(ifig) ],options_);
            else
                if options_.opt_gsa.pprior
                    dyn_saveas(hh,[OutDir '/' fname_ '_rmse_prior_lnpost',int2str(ifig)],options_);
                else
                    dyn_saveas(hh,[OutDir '/' fname_ '_rmse_mc_lnpost',int2str(ifig)],options_);
                end
            end
            if ~options_.nodisplay
                close(hh);
            end
        end
    end
    end
    
    param_names='';
    for j=1:npar+nshock,
        param_names=char(param_names, bayestopt_.name{j});
    end
    param_names=param_names(2:end,:);
    
    disp(' ')
    disp('RMSE over the MC sample:')
    disp('            min yr RMSE    max yr RMSE')
    for j=1:size(vvarvecm,1),
        disp([vvarvecm(j,:), sprintf('%15.5g',[(min(rmse_MC(:,j))) [(max(rmse_MC(:,j)))]])])
    end
    invar = find( std(rmse_MC)./mean(rmse_MC)<=0.0001 );
    if ~isempty(invar)
        disp(' ')
        disp(' ')
        disp('RMSE is not varying significantly over the MC sample for the following variables:')
        disp(vvarvecm(invar,:))
        disp('These variables are excluded from SA')
        disp('[Unless you treat these series as exogenous, there is something wrong in your estimation !]')
    end
    ivar = find( std(rmse_MC)./mean(rmse_MC)>0.0001 );
    vvarvecm=vvarvecm(ivar,:);
    rmse_MC=rmse_MC(:,ivar);
    
    disp(' ')
    % if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
    disp(['Sample filtered the ',num2str(pfilt*100),'% best RMSE''s for each observed series ...' ])
    % else
    %   disp(['Sample filtered the best RMSE''s smaller than RMSE at the posterior mean ...' ])
    % end
    % figure, boxplot(rmse_MC)
    % set(gca,'xticklabel',vvarvecm)
    % saveas(gcf,[fname_,'_SA_RMSE'])
    
    disp(' ')
    disp(' ')
    disp('RMSE ranges after filtering:')
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
        disp(['             best ',num2str(pfilt*100),'% filtered             remaining 90%'])
        disp(['             min            max            min            max            posterior mode'])
    else
        disp(['             best  filtered             remaining '])
        disp(['             min            max            min            max            posterior mean'])
    end
    for j=1:size(vvarvecm,1),
        disp([vvarvecm(j,:), sprintf('%15.5g',[min(rmse_MC(ixx(1:nfilt0(j),j),j)) ...
            max(rmse_MC(ixx(1:nfilt0(j),j),j))  ...
            min(rmse_MC(ixx(nfilt0(j)+1:end,j),j)) ...
            max(rmse_MC(ixx(nfilt0(j)+1:end,j),j)) ...
            rmse_txt(j)])])
        %   disp([vvarvecm(j,:), sprintf('%15.5g',[min(logpo2(ixx(1:nfilt,j))) ...
        %         max(logpo2(ixx(1:nfilt,j)))  ...
        %         min(logpo2(ixx(nfilt+1:end,j))) ...
        %         max(logpo2(ixx(nfilt+1:end,j)))])])
    end
    
    SP=zeros(npar+nshock,size(vvarvecm,1));
    for j=1:size(vvarvecm,1),
        ns=find(PP(:,j)<alpha);
        SP(ns,j)=ones(size(ns));
        SS(:,j)=SS(:,j).*SP(:,j);
    end
    
    for j=1:npar+nshock, %estim_params_.np,
        nsp(j)=length(find(SP(j,:)));
    end
    snam0=param_names(find(nsp==0),:);
    snam1=param_names(find(nsp==1),:);
    snam2=param_names(find(nsp>1),:);
    snam=param_names(find(nsp>0),:);
    % snam0=bayestopt_.name(find(nsp==0));
    % snam1=bayestopt_.name(find(nsp==1));
    % snam2=bayestopt_.name(find(nsp>1));
    % snam=bayestopt_.name(find(nsp>0));
    nsnam=(find(nsp>1));
    
    disp(' ')
    disp(' ')
    disp('These parameters do not affect significantly the fit of ANY observed series:')
    disp(snam0)
    disp(' ')
    disp('These parameters affect ONE single observed series:')
    disp(snam1)
    disp(' ')
    disp('These parameters affect MORE THAN ONE observed series: trade off exists!')
    disp(snam2)
    
    
    %pnam=bayestopt_.name(end-estim_params_.np+1:end);
    pnam=bayestopt_.name;
    
    % plot trade-offs
    if ~options_.nograph
    a00=jet(size(vvarvecm,1));
    for ix=1:ceil(length(nsnam)/5),
        hh = dyn_figure(options_);
        for j=1+5*(ix-1):min(size(snam2,1),5*ix),
            subplot(2,3,j-5*(ix-1))
            %h0=cumplot(x(:,nsnam(j)+nshock));
            h0=cumplot(x(:,nsnam(j)));
            set(h0,'color',[0 0 0])
            hold on,
            npx=find(SP(nsnam(j),:)==0);
            %a0=jet(nsp(nsnam(j)));
%             a0=a00(np,:);
            for i=1:size(vvarvecm,1),
                %h0=cumplot(x(ixx(1:nfilt,np(i)),nsnam(j)+nshock));
%                 h0=cumplot(x(ixx(1:nfilt0(np(i)),np(i)),nsnam(j)));
                if any(npx==i),
                    h0=cumplot(x(ixx(1:nfilt0(i),i),nsnam(j))*NaN);
                else
                    h0=cumplot(x(ixx(1:nfilt0(i),i),nsnam(j)));
            end
                set(h0,'color',a00(i,:))
            end
            ydum=get(gca,'ylim');
            %xdum=xparam1(nshock+nsnam(j));
            if exist('xparam1')
                xdum=xparam1(nsnam(j));
                h1=plot([xdum xdum],ydum);
                set(h1,'color',[0.85 0.85 0.85],'linewidth',2)
            end
            xlabel('')
            title([pnam{nsnam(j)}],'interpreter','none')
        end
        %subplot(3,2,6)
        if exist('OCTAVE_VERSION'),
            legend(char('base',vvarvecm),'location','eastoutside');
        else
            h0=legend(char('base',vvarvecm),0);
            set(h0,'fontsize',6,'position',[0.7 0.1 0.2 0.3],'interpreter','none');
        end
        %h0=legend({'base',vnam{np}}',0);
        %set(findobj(get(h0,'children'),'type','text'),'interpreter','none')
        if options_.opt_gsa.ppost
            dyn_saveas(hh,[ OutDir '/' fname_ '_rmse_post_' int2str(ix)],options_);
        else
            if options_.opt_gsa.pprior
                dyn_saveas(hh,[OutDir '/' fname_ '_rmse_prior_' int2str(ix) ],options_);
            else
                dyn_saveas(hh,[OutDir '/' fname_ '_rmse_mc_' int2str(ix)],options_);
            end
        end
    end
    close all
    end
    
    for j=1:size(SP,2),
        nsx(j)=length(find(SP(:,j)));
    end
    
    disp(' ')
    disp(' ')
    disp('Sensitivity table (significance and direction):')
    vav=char(zeros(1, size(param_names,2)+3 ));
    ibl = 12-size(vvarvecm,2);
    for j=1:size(vvarvecm,1),
        vav = [vav, char(zeros(1,ibl)),vvarvecm(j,:)];
    end
    disp(vav)
    for j=1:npar+nshock, %estim_params_.np,
        %disp([param_names(j,:), sprintf('%8.5g',SP(j,:))])
        disp([param_names(j,:),'   ', sprintf('%12.3g',PP(j,:))])
        disp([char(zeros(1, size(param_names,2)+3 )),sprintf('    (%6g)',SS(j,:))])
    end
    
    
    disp(' ')
    disp(' ')
    disp('Starting bivariate analysis:')
    
    for i=1:size(vvarvecm,1)
        if options_.opt_gsa.ppost
            fnam = ['rmse_post_',deblank(vvarvecm(i,:))];
        else
            if options_.opt_gsa.pprior
                fnam = ['rmse_prior_',deblank(vvarvecm(i,:))];
            else
                fnam = ['rmse_mc_',deblank(vvarvecm(i,:))];
            end
        end
        stab_map_2(x(ixx(1:nfilt0(i),i),:),alpha2,pvalue,fnam, OutDir);
        
        %     [pc,latent,explained] = pcacov(c0);
        %     %figure, bar([explained cumsum(explained)])
        %     ifig=0;
        %     j2=0;
        %     for j=1:npar+nshock,
        %         i2=find(abs(pc(:,j))>alphaPC);
        %         if ~isempty(i2),
        %             j2=j2+1;
        %             if mod(j2,12)==1,
        %                 ifig=ifig+1;
        %                 figure('name',['PCA of the filtered sample ',deblank(vvarvecm(i,:)),' ',num2str(ifig)]),
        %             end
        %             subplot(3,4,j2-(ifig-1)*12)
        %             bar(pc(i2,j)),
        %             set(gca,'xticklabel',bayestopt_.name(i2)),
        %             set(gca,'xtick',[1:length(i2)])
        %             title(['PC ',num2str(j),'. Explained ',num2str(explained(j)),'%'])
        %         end
        %         if (mod(j2,12)==0 | j==(npar+nshock)) & j2,
        %             saveas(gcf,[fname_,'_SA_PCA_',deblank(vvarvecm(i,:)),'_',int2str(ifig)])
        %         end
        %     end
        %     close all
    end
    
end