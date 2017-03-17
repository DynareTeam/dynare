function [rmse_MC, ixx] = filt_mc_(OutDir,options_gsa_,dataset_,dataset_info)
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
% marco.ratto@ec.europa.eu 

% Copyright (C) 2012-2016 European Commission
% Copyright (C) 2012-2016 Dynare Team
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
if options_.TeX
    vvarvecm_tex = options_gsa_.var_rmse_tex;
else
    vvarvecm_tex = [];
end
loadSA   = options_gsa_.load_rmse;
pfilt    = options_gsa_.pfilt_rmse;
alpha    = options_gsa_.alpha_rmse;
% alpha2   = options_gsa_.alpha2_rmse;
alpha2 = 0;
pvalue   = options_gsa_.alpha2_rmse;
istart   = max(2,options_gsa_.istart_rmse);
alphaPC  = 0.5;

fname_ = M_.fname;
lgy_ = M_.endo_names;
dr_ = oo_.dr;

skipline(2)
disp('Starting sensitivity analysis')
disp('for the fit of EACH observed series ...')
skipline()
if ~options_.nograph,
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
        delete([OutDir,filesep,a(j).name])
    end,
    if strmatch([fname_,tmp1],a(j).name),
        disp(a(j).name)
        delete([OutDir,filesep,a(j).name])
    end,
end
disp('done !')
end

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
if loadSA,
    tmplist =load([OutDir,filesep,fnamtmp, '.mat'],'vvarvecm');
    if isempty(fieldnames(tmplist)),
        disp('WARNING: cannot load results since the list of variables used is not present in the mat file')
        loadSA=0;
    elseif ~isequal(tmplist.vvarvecm,vvarvecm)
        disp('WARNING: cannot load results since the list of variables in the mat file differs from the one requested.')
        loadSA=0;
    end
end
if ~loadSA,
    if exist('xparam1','var')
        M_ = set_all_parameters(xparam1,estim_params_,M_);
        ys_mode=steady_(M_,options_,oo_);
    end
    if exist('xparam1_mean','var')
        M_ = set_all_parameters(xparam1_mean,estim_params_,M_);
        ys_mean=steady_(M_,options_,oo_);
    end
    %   eval(options_.datafile)
    Y = transpose(dataset_.data);
    gend = dataset_.nobs;
    data_index = dataset_info.missing.aindex;
    missing_value = dataset_info.missing.state;
    for jx=1:gend
        data_indx(jx,data_index{jx})=true; 
    end
    %stock_gend=data_info.gend;
    %stock_data = data_info.data;
    load([DirectoryName filesep M_.fname '_data.mat']);
    filfilt = dir([DirectoryName filesep M_.fname '_filter_step_ahead*.mat']);
    temp_smooth_file_list = dir([DirectoryName filesep M_.fname '_smooth*.mat']);
    jfile=0;
    for j=1:length(temp_smooth_file_list),
        if isempty(strfind(temp_smooth_file_list(j).name,'smoothed')),
            jfile=jfile+1;
            filsmooth(jfile)=temp_smooth_file_list(j);
        end
    end
    filupdate = dir([DirectoryName filesep M_.fname '_update*.mat']);
    filparam = dir([DirectoryName filesep M_.fname '_param*.mat']);
    x=[];
    logpo2=[];
    sto_ys=[];
    for j=1:length(filparam),
        %load([DirectoryName filesep M_.fname '_param',int2str(j),'.mat']);
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
        skipline()
        disp('Computing RMSE''s...')
        for i=1:size(vvarvecm,1),
            vj=deblank(vvarvecm(i,:));
            
            jxj(i) = strmatch(vj,lgy_(dr_.order_var,:),'exact');
            js(i) = strmatch(vj,lgy_,'exact');
            yss(i,:,:)=repmat(sto_ys(:,js(i))',[gend,1]);
        end
        if exist('xparam1','var')
            [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
            y0 = reshape( squeeze(aK(1,jxj,1:gend)),[gend length(jxj)]);% + kron(ys_mode(js),ones(1,gend)));
            yobs = transpose( ahat(jxj,:));% + kron(ys_mode(js),ones(1,gend)));
            rmse_mode = sqrt(mean((yobs(istart:end,:)-y0(istart:end,:)).^2));
            r2_mode = 1-sum((yobs(istart:end,:)-y0(istart:end,:)).^2)./sum(yobs(istart:end,:).^2);
        end
        y0=-yss;
        nbb=0;
        for j=1:length(filfilt),
            load([DirectoryName filesep M_.fname '_filter_step_ahead',num2str(j),'.mat']);
            nb = size(stock,4);
            y0(:,:,nbb+1:nbb+nb)=y0(:,:,nbb+1:nbb+nb)+reshape(stock(1,js,1:gend,:),[length(js) gend nb]);
            nbb=nbb+nb;
            clear stock;
        end
        yobs=-yss;
        nbb=0;
        for j=1:length(filupdate),
            load([DirectoryName filesep M_.fname '_update',num2str(j),'.mat']);
            nb = size(stock,3);
            yobs(:,:,nbb+1:nbb+nb)=yobs(:,:,nbb+1:nbb+nb)+reshape(stock(js,1:gend,:),[length(js) gend nb]);
            nbb=nbb+nb;
            clear stock;
        end
        y0M=mean(y0,2);
        rmse_MC=zeros(nruns,length(js));
        r2_MC=zeros(nruns,length(js));
        for j=1:nruns,
            rmse_MC(j,:) = sqrt(mean((yobs(:,istart:end,j)'-y0(:,istart:end,j)').^2));
            r2_MC(j,:) = 1-mean((yobs(:,istart:end,j)'-y0(:,istart:end,j)').^2)./mean((yobs(:,istart:end,j)').^2);
        end
        if exist('xparam1_mean','var')
            %eval(['rmse_pmean(i) = sqrt(mean((',vj,'(fobs-1+istart:fobs-1+nobs)-y0M(istart:end-1)).^2));'])
            [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(xparam1_mean,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
            y0 = reshape( squeeze(aK(1,jxj,1:gend)),[gend length(jxj)]);% + kron(ys_mean(js),ones(1,gend)));
            yobs = transpose( ahat(jxj,:));% + kron(ys_mean(js),ones(1,gend)));
            rmse_pmean = sqrt(mean((yobs(istart:end,:)-y0(istart:end,:)).^2));
            r2_pmean = 1-mean((yobs(istart:end,:)-y0(istart:end,:)).^2)./mean(yobs(istart:end,:).^2);
        end
        clear stock_filter;
    end
    for j=1:nruns,
        lnprior(j,1) = priordens(x(j,:)',bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    end
    likelihood=logpo2(:)-lnprior(:);
    disp('... done!')
    
    if options_.opt_gsa.ppost
        save([OutDir,filesep,fnamtmp,'.mat'], 'x', 'logpo2', 'likelihood', 'rmse_MC', 'r2_MC', 'vvarvecm')
        if exist('xparam1_mean','var')
            save([OutDir,filesep,fnamtmp, '.mat'], 'rmse_pmean', 'r2_pmean','-append')
        end
        if exist('xparam1','var')
            save([OutDir,filesep,fnamtmp,'.mat'], 'rmse_mode', 'r2_mode','-append')
        end
    else
        if options_.opt_gsa.lik_only
            save([OutDir,filesep,fnamtmp, '.mat'], 'x', 'logpo2','likelihood', '-append')
        else
            save([OutDir,filesep,fnamtmp, '.mat'], 'x', 'logpo2','likelihood', 'rmse_MC', 'r2_MC', 'vvarvecm','-append')
            if exist('xparam1_mean','var')
                save([OutDir,filesep,fnamtmp, '.mat'], 'rmse_pmean', 'r2_pmean','-append')
            end
            if exist('xparam1','var')
                save([OutDir,filesep,fnamtmp,'.mat'], 'rmse_mode', 'r2_mode','-append')
            end
        end
    end
else
    if options_.opt_gsa.lik_only && options_.opt_gsa.ppost==0
        load([OutDir,filesep,fnamtmp, '.mat'],'x','logpo2','likelihood');
    else
        load([OutDir,filesep,fnamtmp, '.mat'],'x','logpo2','likelihood','rmse_MC','rmse_mode','rmse_pmean', 'r2_MC', 'vvarvecm', 'r2_mode','r2_pmean');
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

%% visual scatter analysis!
if options_.opt_gsa.ppost
    tmp_title='R2 Posterior:';
    atitle='R2 Posterior:';
    asname='r2_post';
else
    if options_.opt_gsa.pprior
        tmp_title='R2 Prior:';
        atitle='R2 Prior:';
        asname='r2_prior';
    else
        tmp_title='R2 MC:';
        atitle='R2 MC:';
        asname='r2_mc';
    end
end
options_scatter.param_names = vvarvecm;
options_scatter.param_names_tex = vvarvecm_tex;
options_scatter.fname_ = fname_;
options_scatter.OutputDirectoryName = OutDir;
options_scatter.amcf_name = asname;
options_scatter.amcf_title = atitle;
options_scatter.title = tmp_title;
scatter_analysis(r2_MC, x,options_scatter, options_);
%% end of visual scatter analysis


if ~options_.opt_gsa.ppost && options_.opt_gsa.lik_only
    if options_.opt_gsa.pprior
        anam='rmse_prior_post';
        atitle='RMSE prior: Log Posterior Kernel';
    else
        anam='rmse_mc_post';
        atitle='RMSE MC: Log Posterior Kernel';
    end
    options_mcf.pvalue_ks = alpha;
    options_mcf.pvalue_corr = pvalue;
    options_mcf.alpha2 = alpha2;
    if options_.TeX
        [pnames,pnames_tex]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        options_mcf.param_names = char(pnames);
        options_mcf.param_names_tex = char(pnames_tex);
    else
        [pnames]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        options_mcf.param_names = char(pnames);
        options_mcf.param_names_tex = [];
    end
    options_mcf.fname_ = fname_;
    options_mcf.OutputDirectoryName = OutDir;
    options_mcf.amcf_name = anam;
    options_mcf.amcf_title = atitle;
    options_mcf.title = atitle;
    options_mcf.beha_title = 'better posterior kernel';
    options_mcf.nobeha_title = 'worse posterior kernel';
    mcf_analysis(x, ipost(1:nfilt), ipost(nfilt+1:end), options_mcf, options_);
    
    if options_.opt_gsa.pprior
        anam = 'rmse_prior_lik';
        atitle = 'RMSE prior: Log Likelihood Kernel';
    else
        anam='rmse_mc_lik';
        atitle = 'RMSE MC: Log Likelihood Kernel';
    end
    options_mcf.amcf_name = anam;
    options_mcf.amcf_title = atitle;
    options_mcf.title = atitle;
    options_mcf.beha_title = 'better likelihood';
    options_mcf.nobeha_title = 'worse likelihood';
    mcf_analysis(x, ilik(1:nfilt), ilik(nfilt+1:end), options_mcf, options_);
    
else
    if options_.opt_gsa.ppost,
        rmse_txt=rmse_pmean;
        r2_txt=r2_pmean;
    else
        if options_.opt_gsa.pprior || ~exist('rmse_pmean'),
            if exist('rmse_mode'),
                rmse_txt=rmse_mode;
                r2_txt=r2_mode;
            else
                rmse_txt=NaN(1,size(rmse_MC,2));
                r2_txt=NaN(1,size(r2_MC,2));
            end
        else
            %nfilt0(i)=length(find(rmse_MC(:,i)<rmse_pmean(i)));
            rmse_txt=rmse_pmean;
            r2_txt=r2_pmean;
        end
    end
    for i=1:size(vvarvecm,1),
        [dum, ixx(:,i)]=sort(rmse_MC(:,i));
    end
    PP=ones(npar+nshock,size(vvarvecm,1));
    PPV=ones(size(vvarvecm,1),size(vvarvecm,1),npar+nshock);
    SS=zeros(npar+nshock,size(vvarvecm,1));
    for j=1:npar+nshock,
        for i=1:size(vvarvecm,1),
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
        for i=1:size(vvarvecm,1),
            for l=1:size(vvarvecm,1),
                if l~=i && PP(j,i)<alpha && PP(j,l)<alpha
                    [H,P,KSSTAT] = smirnov(x(ixx(1:nfilt0(i),i),j),x(ixx(1:nfilt0(l),l),j), alpha);
                    %[H1,P1,KSSTAT1] = smirnov(x(ixx(1:nfilt0(i),i),j),x(:,j), alpha);
                    %                 PP(j,i)=min(P,PP(j,i));
                    %                 PP(j,l)=min(P,PP(j,l));
                    %if P<P1
                    %                     if SS(j,i)*SS(j,l)
                    PPV(i,l,j) = P;
                    %                     end
                elseif l==i
                    PPV(i,l,j) = PP(j,i);
                end
            end
        end
    end
    if ~options_.nograph,
        ifig=0;
        for i=1:size(vvarvecm,1),
            if options_.opt_gsa.ppost
                temp_name='RMSE Posterior: Log Prior';
            else
                if options_.opt_gsa.pprior
                    temp_name='RMSE Prior: Log Prior';
                else
                    temp_name='RMSE MC: Log Prior';
                end
            end
            if mod(i,9)==1,
                ifig=ifig+1;
                hh=dyn_figure(options_,'name',[temp_name,' ',int2str(ifig)]);
            end
            subplot(3,3,i-9*(ifig-1))
            h=cumplot(lnprior(ixx(1:nfilt0(i),i)));
            set(h,'color','blue','linewidth',2)
            hold on, h=cumplot(lnprior);
            set(h,'color','k','linewidth',1)
            h=cumplot(lnprior(ixx(nfilt0(i)+1:end,i)));
            set(h,'color','red','linewidth',2)
            title(vvarvecm(i,:),'interpreter','none')
            if mod(i,9)==0 || i==size(vvarvecm,1)
                if ~isoctave
                    annotation('textbox', [0.1,0,0.35,0.05],'String', 'Log-prior for BETTER R2','Color','Blue','horizontalalignment','center');
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Log-prior for WORSE R2', 'Color','Red','horizontalalignment','center');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh,[OutDir filesep fname_ '_rmse_post_lnprior',int2str(ifig)],options_);
                    if options_.TeX
                    create_TeX_loader(options_,[OutDir '/' fname_ '_rmse_post_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_post_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1))
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_prior_lnprior',int2str(ifig) ],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir '/' fname_ '_rmse_prior_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_prior_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1))
                    end
                    else
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_mc_lnprior',int2str(ifig) ],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir '/' fname_ '_rmse_mc_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_mc_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1))
                    end
                    end
                end
            end
        end
        ifig=0;
        for i=1:size(vvarvecm,1),
            if options_.opt_gsa.ppost
                temp_name='RMSE Posterior: Log Likelihood';
            else
                if options_.opt_gsa.pprior
                    temp_name='RMSE Prior: Log Likelihood';
                else
                    temp_name='RMSE MC: Log Likelihood';
                end
            end
            if mod(i,9)==1,
                ifig=ifig+1;
                hh = dyn_figure(options_,'Name',[temp_name,' ',int2str(ifig)]);
            end
            subplot(3,3,i-9*(ifig-1))
            h=cumplot(likelihood(ixx(1:nfilt0(i),i)));
            set(h,'color','blue','linewidth',2)
            hold on, h=cumplot(likelihood);
            set(h,'color','k','linewidth',1)
            h=cumplot(likelihood(ixx(nfilt0(i)+1:end,i)));
            set(h,'color','red','linewidth',2)
            title(vvarvecm(i,:),'interpreter','none')
            if options_.opt_gsa.ppost==0,
                set(gca,'xlim',[min( likelihood(ixx(1:nfilt0(i),i)) ) max( likelihood(ixx(1:nfilt0(i),i)) )])
            end
            if mod(i,9)==0 || i==size(vvarvecm,1)
                if ~isoctave
                    annotation('textbox', [0.1,0,0.35,0.05],'String', 'Log-likelihood for BETTER R2','Color','Blue','horizontalalignment','center');
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Log-likelihood for WORSE R2', 'Color','Red','horizontalalignment','center');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh,[OutDir filesep fname_ '_rmse_post_lnlik',int2str(ifig) ],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_post_lnprior',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_post_lnprior',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_prior_lnlik',int2str(ifig)],options_);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_lnlik',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_prior_lnlik',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    else
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_mc_lnlik',int2str(ifig) ],options_);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_lnlik',int2str(ifig) ],ifig,[temp_name,' ',int2str(ifig)],'rmse_mc_lnlik',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    end
                end
            end
        end
        ifig=0;
        for i=1:size(vvarvecm,1),
            if options_.opt_gsa.ppost
                temp_name='RMSE Posterior: Log Posterior';
            else
                if options_.opt_gsa.pprior
                    temp_name='RMSE Prior: Log Posterior';
                else
                    temp_name='RMSE MC: Log Posterior';
                end
            end
            if mod(i,9)==1,
                ifig=ifig+1;
                hh = dyn_figure(options_,'Name',[temp_name,' ',int2str(ifig)]);
            end
            subplot(3,3,i-9*(ifig-1))
            h=cumplot(logpo2(ixx(1:nfilt0(i),i)));
            set(h,'color','blue','linewidth',2)
            hold on, h=cumplot(logpo2);
            set(h,'color','k','linewidth',1)
            h=cumplot(logpo2(ixx(nfilt0(i)+1:end,i)));
            set(h,'color','red','linewidth',2)
            title(vvarvecm(i,:),'interpreter','none')
            if options_.opt_gsa.ppost==0,
                set(gca,'xlim',[min( logpo2(ixx(1:nfilt0(i),i)) ) max( logpo2(ixx(1:nfilt0(i),i)) )])
            end
            if mod(i,9)==0 || i==size(vvarvecm,1)
                if ~isoctave
                    annotation('textbox', [0.1,0,0.35,0.05],'String', 'Log-posterior for BETTER R2','Color','Blue','horizontalalignment','center');
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Log-posterior for WORSE R2', 'Color','Red','horizontalalignment','center');
                end
                if options_.opt_gsa.ppost
                    dyn_saveas(hh,[OutDir filesep fname_ '_rmse_post_lnpost',int2str(ifig) ],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_post_lnpost',int2str(ifig) ],ifig,[temp_name,' ',int2str(ifig)],'rmse_post_lnpost',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_prior_lnpost',int2str(ifig)],options_);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_lnpost',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_prior_lnpost',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    else
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_mc_lnpost',int2str(ifig)],options_);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_lnpost',int2str(ifig)],ifig,[temp_name,' ',int2str(ifig)],'rmse_mc_lnpost',options_.figures.textwidth*min((i-9*(ifig-1))/3,1));
                        end
                    end
                end
            end
        end
    end
    
    if options_.TeX
        [pnames,pnames_tex]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        param_names = char(pnames);
        param_names_tex = char(pnames_tex);
    else
        [pnames]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_);
        param_names = char(pnames);
        param_names_tex = [];
    end

    
    skipline()
    title_string='RMSE over the MC sample:';
    data_mat=[min(rmse_MC)' max(rmse_MC)'];
    headers=strvcat('Variable','min yr RMSE','max yr RMSE');
    dyntable(options_,title_string,headers,vvarvecm,data_mat, 0, 15, 5);
    if options_.TeX
        headers_tex=strvcat('\text{Variable}','\text{min yr RMSE}','\text{max yr RMSE}');
        dyn_latex_table(M_,options_,title_string,'RMSE_MC',headers_tex,vvarvecm_tex,data_mat,0,15,5);
    end
    
    invar = find( std(rmse_MC)./mean(rmse_MC)<=0.0001 );
    if ~isempty(invar)
        skipline(2)
        disp('RMSE is not varying significantly over the MC sample for the following variables:')
        disp(vvarvecm(invar,:))
        disp('These variables are excluded from SA')
        disp('[Unless you treat these series as exogenous, there is something wrong in your estimation !]')
    end
    vvarvecm0=vvarvecm;

    ivar = find( std(rmse_MC)./mean(rmse_MC)>0.0001 );
    vvarvecm=vvarvecm(ivar,:);
    rmse_MC=rmse_MC(:,ivar);
    
    skipline()
    % if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
    disp(['Sample filtered the ',num2str(pfilt*100),'% best RMSE''s for each observed series ...' ])
    % else
    %   disp(['Sample filtered the best RMSE''s smaller than RMSE at the posterior mean ...' ])
    % end
    % figure, boxplot(rmse_MC)
    % set(gca,'xticklabel',vvarvecm)
    % saveas(gcf,[fname_,'_SA_RMSE'])
    
    skipline(2)
    disp('RMSE ranges after filtering:')
    title_string='RMSE ranges after filtering:';
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
        headers=strvcat('Variable','min','max','min','max','posterior mode');
        headers_tex=strvcat('\text{Variable}','\text{min}','\text{max}','\text{min}','\text{max}','\text{posterior mode}');
    else
        headers=strvcat('Variable','min','max','min','max','posterior mean');
        headers_tex=strvcat('\text{Variable}','\text{min}','\text{max}','\text{min}','\text{max}','\text{posterior mean}');
    end
    data_mat=NaN(size(vvarvecm,1),5);
    for j=1:size(vvarvecm,1),
        data_mat(j,:)=[min(rmse_MC(ixx(1:nfilt0(j),j),j)) ...
            max(rmse_MC(ixx(1:nfilt0(j),j),j))  ...
            min(rmse_MC(ixx(nfilt0(j)+1:end,j),j)) ...
            max(rmse_MC(ixx(nfilt0(j)+1:end,j),j)) ...
            rmse_txt(j)];
    end
    %get formatting for additional header line
    val_width=15;
    val_precis=5;
    label_width = max(size(deblank(char(headers(1,:),vvarvecm)),2)+2,0);
    label_format_leftbound  = sprintf('%%-%ds',label_width);
    if all(~isfinite(data_mat))
        values_length = 4;
    else
        values_length = max(ceil(max(max(log10(abs(data_mat(isfinite(data_mat))))))),1)+val_precis+1;
    end
    if any(data_mat) < 0 %add one character for minus sign
        values_length = values_length+1;
    end
    
    headers_length = max(size(deblank(headers(2:end,:)),2));
    if ~isempty(val_width)
        val_width = max(max(headers_length,values_length)+2,val_width);
    else
        val_width = max(headers_length,values_length)+2;
    end
    value_format  = sprintf('%%%d.%df',val_width,val_precis);
    header_string_format  = sprintf('%%%ds',val_width);
    
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
        optional_header=sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','',['best ',num2str(pfilt*100),'% filtered'],'','remaining 90%');
    else
        optional_header=sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','','best  filtered','','remaining');
    end
    dyntable(options_,title_string,headers,vvarvecm,data_mat, 0, val_width, val_precis,optional_header);
    if options_.TeX
        if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
            optional_header={[' & \multicolumn{2}{c}{best ',num2str(pfilt*100),' filtered} & \multicolumn{2}{c}{remaining 90\%}\\']};
        else
            optional_header={[' & \multicolumn{2}{c}{best filtered} & \multicolumn{2}{c}{remaining}\\']};
        end
        dyn_latex_table(M_,options_,title_string,'RMSE_ranges_after_filtering',headers_tex,vvarvecm_tex,data_mat,0,val_width,val_precis,optional_header);
    end
    
    %%%%% R2 table
    vvarvecm=vvarvecm0;
    skipline()
    title_string='R2 over the MC sample:';
    data_mat=[min(r2_MC)' max(r2_MC)'];
    headers=strvcat('Variable','min yr R2','max yr R2');
    dyntable(options_,title_string,headers,vvarvecm,data_mat, 0, 15, 5);
    if options_.TeX
        headers_tex=strvcat('\text{Variable}','\text{min yr R2}','\text{max yr R2}');
        dyn_latex_table(M_,options_,title_string,'R2_MC',headers_tex,vvarvecm_tex,data_mat,0,15,5);
    end
    
    r2_MC=r2_MC(:,ivar);
    vvarvecm=vvarvecm(ivar,:);
    
    skipline()
    disp(['Sample filtered the ',num2str(pfilt*100),'% best R2''s for each observed series ...' ])
    
    skipline()
    disp('R2 ranges after filtering:')
    title_string='R2 ranges after filtering:';
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
        headers=strvcat('Variable','min','max','min','max','posterior mode');
        headers_tex=strvcat('\text{Variable}','\text{min}','\text{max}','\text{min}','\text{max}','\text{posterior mode}');
    else
        headers=strvcat('Variable','min','max','min','max','posterior mean');
        headers_tex=strvcat('\text{Variable}','\text{min}','\text{max}','\text{min}','\text{max}','\text{posterior mean}');
    end
    data_mat=NaN(size(vvarvecm,1),5);
    for j=1:size(vvarvecm,1),
        data_mat(j,:)=[min(r2_MC(ixx(1:nfilt0(j),j),j)) ...
            max(r2_MC(ixx(1:nfilt0(j),j),j))  ...
            min(r2_MC(ixx(nfilt0(j)+1:end,j),j)) ...
            max(r2_MC(ixx(nfilt0(j)+1:end,j),j)) ...
            r2_txt(j)];
    end
    %get formatting for additional header line
    val_width=15;
    val_precis=5;
    label_width = max(size(deblank(char(headers(1,:),vvarvecm)),2)+2,0);
    label_format_leftbound  = sprintf('%%-%ds',label_width);
    if all(~isfinite(data_mat))
        values_length = 4;
    else
        values_length = max(ceil(max(max(log10(abs(data_mat(isfinite(data_mat))))))),1)+val_precis+1;
    end
    if any(data_mat) < 0 %add one character for minus sign
        values_length = values_length+1;
    end
    
    headers_length = max(size(deblank(headers(2:end,:)),2));
    if ~isempty(val_width)
        val_width = max(max(headers_length,values_length)+2,val_width);
    else
        val_width = max(headers_length,values_length)+2;
    end
    value_format  = sprintf('%%%d.%df',val_width,val_precis);
    header_string_format  = sprintf('%%%ds',val_width);
    
    if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
        optional_header=sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','',['best ',num2str(pfilt*100),'% filtered'],'','remaining 90%');
    else
        optional_header=sprintf([label_format_leftbound,header_string_format,header_string_format,header_string_format,header_string_format],'','','best  filtered','','remaining');
    end
    dyntable(options_,title_string,headers,vvarvecm,data_mat, 0, val_width, val_precis,optional_header);
    if options_.TeX
        if options_.opt_gsa.ppost==0 && options_.opt_gsa.pprior,
            optional_header={[' & \multicolumn{2}{c}{best ',num2str(pfilt*100),' filtered} & \multicolumn{2}{c}{remaining 90\%}\\']};
        else
            optional_header={[' & \multicolumn{2}{c}{best filtered} & \multicolumn{2}{c}{remaining}\\']};
        end
        dyn_latex_table(M_,options_,title_string,'R2_ranges_after_filtering',headers_tex,vvarvecm_tex,data_mat,0,val_width,val_precis,optional_header);
    end
    
    
    %%%%  R2 table
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
    
    skipline(2)
    disp('These parameters do not affect significantly the fit of ANY observed series:')
    disp(snam0)
    skipline()
    disp('These parameters affect ONE single observed series:')
    disp(snam1)
    skipline()
    disp('These parameters affect MORE THAN ONE observed series: trade off exists!')
    disp(snam2)
    
    
    %pnam=bayestopt_.name(end-estim_params_.np+1:end);
    pnam=bayestopt_.name;
    
    % plot trade-offs
    if ~options_.nograph
        a00=jet(size(vvarvecm,1));
        if options_.opt_gsa.ppost
            temp_name='RMSE Posterior Tradeoffs:';
            atitle='RMSE Posterior Map:';
            asname='rmse_post';
        else
            if options_.opt_gsa.pprior
                temp_name='RMSE Prior Tradeoffs:';
                atitle='RMSE Prior Map:';
                asname='rmse_prior';
            else
                temp_name='RMSE MC Tradeoffs:';
                atitle='RMSE MC Map:';
                asname='rmse_mc';
            end
        end
        % now I plot by observed variables
        options_mcf.pvalue_ks = alpha;
        options_mcf.pvalue_corr = pvalue;
        options_mcf.alpha2 = alpha2;
        options_mcf.param_names = param_names;
        options_mcf.param_names_tex = param_names_tex;
        options_mcf.fname_ = fname_;
        options_mcf.OutputDirectoryName = OutDir;
        for iy=1:size(vvarvecm,1),
            options_mcf.amcf_name = [asname '_' deblank(vvarvecm(iy,:)) '_map' ];
            options_mcf.amcf_title = [atitle ' ' deblank(vvarvecm(iy,:))];
            options_mcf.beha_title = ['better fit of ' deblank(vvarvecm(iy,:))];
            options_mcf.nobeha_title = ['worse fit of ' deblank(vvarvecm(iy,:))];
            options_mcf.title = ['the fit of ' deblank(vvarvecm(iy,:))];
            mcf_analysis(x, ixx(1:nfilt0(iy),iy), ixx(nfilt0(iy)+1:end,iy), options_mcf, options_);
        end
        for iy=1:size(vvarvecm,1),
            ipar = find(any(squeeze(PPV(iy,:,:))<alpha));
            for ix=1:ceil(length(ipar)/5),
                hh = dyn_figure(options_,'name',[temp_name,' observed variable ',deblank(vvarvecm(iy,:))]);
                for j=1+5*(ix-1):min(length(ipar),5*ix),
                    subplot(2,3,j-5*(ix-1))
                    %h0=cumplot(x(:,nsnam(j)+nshock));
                    h0=cumplot(x(:,ipar(j)));
                    set(h0,'color',[0 0 0])
                    hold on,
                    iobs=find(squeeze(PPV(iy,:,ipar(j)))<alpha);
                    for i=1:size(vvarvecm,1),
                        %h0=cumplot(x(ixx(1:nfilt,np(i)),nsnam(j)+nshock));
                        %                 h0=cumplot(x(ixx(1:nfilt0(np(i)),np(i)),nsnam(j)));
                        if any(iobs==i) || i==iy,
                            h0=cumplot(x(ixx(1:nfilt0(i),i),ipar(j)));
                            if ~isoctave
                                hcmenu = uicontextmenu;
                                uimenu(hcmenu,'Label',deblank(vvarvecm(i,:)));
                                set(h0,'uicontextmenu',hcmenu)
                            end
                        else
                            h0=cumplot(x(ixx(1:nfilt0(i),i),ipar(j))*NaN);
                        end
                        set(h0,'color',a00(i,:),'linewidth',2)
                    end
                    ydum=get(gca,'ylim');
                    %xdum=xparam1(nshock+nsnam(j));
                    if exist('xparam1')
                        xdum=xparam1(ipar(j));
                        h1=plot([xdum xdum],ydum);
                        set(h1,'color',[0.85 0.85 0.85],'linewidth',2)
                    end
                    xlabel('')
                    title([pnam{ipar(j)}],'interpreter','none')
                end
                %subplot(3,2,6)
                if isoctave
                    legend(char('base',vvarvecm),'location','eastoutside');
                else
                    h0=legend(char('base',vvarvecm));
                    set(h0,'fontsize',6,'position',[0.7 0.1 0.2 0.3],'interpreter','none');
                end
                %h0=legend({'base',vnam{np}}',0);
                %set(findobj(get(h0,'children'),'type','text'),'interpreter','none')
                if options_.opt_gsa.ppost
                    dyn_saveas(hh,[ OutDir filesep fname_ '_rmse_post_' deblank(vvarvecm(iy,:)) '_' int2str(ix)],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[ OutDir filesep fname_ '_rmse_post_' deblank(vvarvecm(iy,:)) '_' int2str(ix)],ix,[temp_name,' observed variable $',deblank(vvarvecm_tex(iy,:)) '$'],['rmse_post_' deblank(vvarvecm(iy,:))],1)
                    end
                else
                    if options_.opt_gsa.pprior
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_prior_' deblank(vvarvecm(iy,:)) '_' int2str(ix) ],options_);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_' deblank(vvarvecm(iy,:)) '_' int2str(ix) ],ix,[temp_name,' observed variable $',deblank(vvarvecm_tex(iy,:)) '$'],['rmse_prior_' deblank(vvarvecm(iy,:))],1)
                        end
                    else
                        dyn_saveas(hh,[OutDir filesep fname_ '_rmse_mc_' deblank(vvarvecm(iy,:)) '_' int2str(ix)],options_);
                        if options_.TeX
                            create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_' deblank(vvarvecm(iy,:)) '_' int2str(ix)],ix,[temp_name,' observed variable $',deblank(vvarvecm_tex(iy,:)) '$'],['rmse_mc_' deblank(vvarvecm(iy,:))],1)
                        end
                    end
                end
            end
        end
        
        % now I plot by individual parameters
        for ix=1:ceil(length(nsnam)/5),
            hh = dyn_figure(options_,'name',[temp_name,' estimated params and shocks ',int2str(ix)]);
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
                        if ~isoctave
                            hcmenu = uicontextmenu;
                            uimenu(hcmenu,'Label',deblank(vvarvecm(i,:)));
                            set(h0,'uicontextmenu',hcmenu)
                        end
                    end
                    set(h0,'color',a00(i,:),'linewidth',2)
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
            if isoctave
                legend(char('base',vvarvecm),'location','eastoutside');
            else
                h0=legend(char('base',vvarvecm));
                set(h0,'fontsize',6,'position',[0.7 0.1 0.2 0.3],'interpreter','none');
            end
            %h0=legend({'base',vnam{np}}',0);
            %set(findobj(get(h0,'children'),'type','text'),'interpreter','none')
            if options_.opt_gsa.ppost
                dyn_saveas(hh,[ OutDir filesep fname_ '_rmse_post_params_' int2str(ix)],options_);
                if options_.TeX
                    create_TeX_loader(options_,[ OutDir filesep fname_ '_rmse_post_params_' int2str(ix)],ix,[temp_name,' estimated params and shocks ',int2str(ix)],'rmse_post_params',1)
                end
            else
                if options_.opt_gsa.pprior
                    dyn_saveas(hh,[OutDir filesep fname_ '_rmse_prior_params_' int2str(ix) ],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_prior_params_' int2str(ix) ],ix,[temp_name,' estimated params and shocks ',int2str(ix)],'rmse_prior_params',1)
                    end
                else
                    dyn_saveas(hh,[OutDir filesep fname_ '_rmse_mc_params_' int2str(ix)],options_);
                    if options_.TeX
                        create_TeX_loader(options_,[OutDir filesep fname_ '_rmse_mc_params_' int2str(ix)],ix,[temp_name,' estimated params and shocks ',int2str(ix)],'rmse_mc_params',1)
                    end
                end
            end
        end
    end
    
%     for j=1:size(SP,2),
%         nsx(j)=length(find(SP(:,j)));
%     end
    
%     skipline(2)
%     disp('Sensitivity table (significance and direction):')
%     vav=char(zeros(1, size(param_names,2)+3 ));
%     ibl = 12-size(vvarvecm,2);
%     for j=1:size(vvarvecm,1),
%         vav = [vav, char(zeros(1,ibl)),vvarvecm(j,:)];
%     end
%     disp(vav)
%     for j=1:npar+nshock, %estim_params_.np,
%         %disp([param_names(j,:), sprintf('%8.5g',SP(j,:))])
%         disp([param_names(j,:),'   ', sprintf('%12.3g',PP(j,:))])
%         disp([char(zeros(1, size(param_names,2)+3 )),sprintf('    (%6g)',SS(j,:))])
%     end
    
    
%     skipline()
%     disp('Starting bivariate analysis:')
%     
%     for i=1:size(vvarvecm,1)
%         if options_.opt_gsa.ppost
%             fnam = ['rmse_post_',deblank(vvarvecm(i,:))];
%         else
%             if options_.opt_gsa.pprior
%                 fnam = ['rmse_prior_',deblank(vvarvecm(i,:))];
%             else
%                 fnam = ['rmse_mc_',deblank(vvarvecm(i,:))];
%             end
%         end
%         stab_map_2(x(ixx(1:nfilt0(i),i),:),alpha2,pvalue,fnam, OutDir,[],[temp_name ' observed variable ' deblank(vvarvecm(i,:))]);
%         
%         %     [pc,latent,explained] = pcacov(c0);
%         %     %figure, bar([explained cumsum(explained)])
%         %     ifig=0;
%         %     j2=0;
%         %     for j=1:npar+nshock,
%         %         i2=find(abs(pc(:,j))>alphaPC);
%         %         if ~isempty(i2),
%         %             j2=j2+1;
%         %             if mod(j2,12)==1,
%         %                 ifig=ifig+1;
%         %                 figure('name',['PCA of the filtered sample ',deblank(vvarvecm(i,:)),' ',num2str(ifig)]),
%         %             end
%         %             subplot(3,4,j2-(ifig-1)*12)
%         %             bar(pc(i2,j)),
%         %             set(gca,'xticklabel',bayestopt_.name(i2)),
%         %             set(gca,'xtick',[1:length(i2)])
%         %             title(['PC ',num2str(j),'. Explained ',num2str(explained(j)),'%'])
%         %         end
%         %         if (mod(j2,12)==0 | j==(npar+nshock)) & j2,
%         %             saveas(gcf,[fname_,'_SA_PCA_',deblank(vvarvecm(i,:)),'_',int2str(ifig)])
%         %         end
%         %     end
%         %     close all
%     end
    
end

function []=create_TeX_loader(options_,figpath,label_number,caption,label_name,scale_factor)
if nargin<6
    scale_factor=1;
end
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([figpath '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by filt_mc_.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',scale_factor,strrep(figpath,'\','/'));
    fprintf(fidTeX,'\\caption{%s.}',caption);
    fprintf(fidTeX,'\\label{Fig:%s:%u}\n',label_name,label_number);
    fprintf(fidTeX,'\\end{figure}\n\n');
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end

function [pnames,pnames_tex]=get_LaTeX_parameter_names(M_,options_,estim_params_,bayestopt_)
np=size(bayestopt_.name,1);
pnames=cell(np,1);
pnames_tex=cell(np,1);    
for ii=1:length(bayestopt_.name)
    if options_.TeX
        [param_name_temp, param_name_tex_temp]= get_the_name(ii,options_.TeX,M_,estim_params_,options_);
        pnames_tex{ii,1} = strrep(param_name_tex_temp,'$','');
        pnames{ii,1} = param_name_temp;
    else
        param_name_temp = get_the_name(ii,options_.TeX,M_,estim_params_,options_);
        pnames{ii,1} = param_name_temp;
    end
end

