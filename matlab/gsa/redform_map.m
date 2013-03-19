function redform_map(dirname,options_gsa_)
%function redform_map(dirname)
% inputs (from opt_gsa structure
% anamendo    = options_gsa_.namendo;
% anamlagendo = options_gsa_.namlagendo;
% anamexo     = options_gsa_.namexo;
% iload       = options_gsa_.load_redform;
% pprior      = options_gsa_.pprior;
% ilog        = options_gsa_.logtrans_redform;
% threshold   = options_gsa_.threshold_redform;
% ksstat      = options_gsa_.ksstat_redform;
% alpha2      = options_gsa_.alpha2_redform;
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


global M_ oo_ estim_params_ options_ bayestopt_

% options_gsa_ = options_.opt_gsa;

anamendo = options_gsa_.namendo;
anamlagendo = options_gsa_.namlagendo;
anamexo = options_gsa_.namexo;
iload = options_gsa_.load_redform;
pprior = options_gsa_.pprior;
ilog = options_gsa_.logtrans_redform;
threshold = options_gsa_.threshold_redform;
% ksstat = options_gsa_.ksstat_redform;
alpha2 = options_gsa_.alpha2_redform;
alpha2=0;
pvalue_ks = options_gsa_.ksstat_redform;
pvalue_corr = options_gsa_.alpha2_redform;

pnames = M_.param_names(estim_params_.param_vals(:,1),:);
if nargin==0,
    dirname='';
end

if pprior
    load([dirname,'/',M_.fname,'_prior'],'lpmat', 'lpmat0', 'istable','T');
    adir=[dirname '/redform_stab'];
else
    load([dirname,'/',M_.fname,'_mc'],'lpmat', 'lpmat0', 'istable','T');
    adir=[dirname '/redform_mc'];
end
if ~exist('T')
    stab_map_(dirname,options_gsa_);
    if pprior
        load([dirname,'/',M_.fname,'_prior'],'T');
    else
        load([dirname,'/',M_.fname,'_mc'],'T');
    end
    if ~exist('T'),
        disp('The model is too large!')
        disp('Reduced form mapping stopped!')
        return
    end
end
if isempty(dir(adir))
    mkdir(adir)
end
adir0=pwd;
%cd(adir)

nspred=size(T,2)-M_.exo_nbr;
x0=lpmat(istable,:);
xx0=lpmat0(istable,:);
[kn, np]=size(x0);
offset = length(bayestopt_.pshape)-np;
if options_gsa_.prior_range,
    pshape=5*(ones(np,1));
    pd =  [NaN(np,1) NaN(np,1) bayestopt_.lb(offset+1:end) bayestopt_.ub(offset+1:end)];
else
    pshape = bayestopt_.pshape(offset+1:end);
    pd =  [bayestopt_.p6(offset+1:end) bayestopt_.p7(offset+1:end) bayestopt_.p3(offset+1:end) bayestopt_.p4(offset+1:end)];
end

nsok = length(find(M_.lead_lag_incidence(M_.maximum_lag,:)));
clear lpmat lpmat0
js=0;
for j=1:size(anamendo,1)
    namendo=deblank(anamendo(j,:));
    iendo=strmatch(namendo,M_.endo_names(oo_.dr.order_var,:),'exact');
    ifig=0;
    iplo=0;
    for jx=1:size(anamexo,1)
        namexo=deblank(anamexo(jx,:));
        iexo=strmatch(namexo,M_.exo_names,'exact');
        disp(' ')
        disp(['[', namendo,' vs. ',namexo,']'])

        
        if ~isempty(iexo),
            %y0=squeeze(T(iendo,iexo+nspred,istable));
            y0=squeeze(T(iendo,iexo+nspred,:));
            if (max(y0)-min(y0))>1.e-10,
                if mod(iplo,9)==0 && isempty(threshold) && ~options_.nograph,
                    ifig=ifig+1;
                    hfig = dyn_figure(options_,'name',[namendo,' vs. shocks ',int2str(ifig)]);
                    iplo=0;
                end
                iplo=iplo+1;
                js=js+1;
                xdir0 = [adir,'/',namendo,'_vs_', namexo];
                if ilog==0,
                    if isempty(threshold)
                        if isempty(dir(xdir0))
                            mkdir(xdir0)
                        end
                        si(:,js) = redform_private(x0, y0, pshape, pd, iload, pnames, namendo, namexo, xdir0, options_gsa_);
                    else
                        iy=find( (y0>threshold(1)) & (y0<threshold(2)));
                        iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
                        xdir = [xdir0,'_cut'];
                        if isempty(dir(xdir))
                            mkdir(xdir)
                        end
                        if ~options_.nograph,
                            hf=dyn_figure(options_); hist(y0,30), title([namendo,' vs. ', namexo])
                            dyn_saveas(hf,[xdir,'/', namendo,'_vs_', namexo],options_);
                            if ~options_.nodisplay
                                close(hf);
                            end
                        end
                        %             if ~isempty(iy),
                        %               si(:,js) = redform_private(x0(iy,:), y0(iy), pshape, pd, iload, pnames, namendo, namexo, xdir, options_gsa_);
                        %             else
                        si(:,js) = NaN(np,1);
                        %             end
                        if ~isempty(iy) && ~isempty(iyc)
                            delete([xdir, '/*cut*.*'])
                            [proba, dproba] = stab_map_1(x0, iy, iyc, 'cut',0);
                            %             indsmirnov = find(dproba>ksstat);
                            indsmirnov = find(proba<pvalue_ks);
                            for jp=1:length(indsmirnov),
                                disp([M_.param_names(estim_params_.param_vals(indsmirnov(jp),1),:),'   d-stat = ', num2str(dproba(indsmirnov(jp)),'%1.3f'),'   p-value = ', num2str(proba(indsmirnov(jp)),'%1.3f')])
                            end
                            disp(' ');
                            stab_map_1(x0, iy, iyc, 'cut',1,indsmirnov,xdir);
                            stab_map_2(x0(iy,:),alpha2,pvalue_corr,'cut',xdir)
                            stab_map_2(x0(iyc,:),alpha2,pvalue_corr,'trim',xdir)
                            lpmat=x0(iy,:);
                            lpmat0=xx0(iy,:);
                            istable=[1:length(iy)];
                            save([xdir,filesep,'threshold.mat'],'lpmat','lpmat0','istable')
                            clear lpmat lpmat0 istable
                        end
                    end
                else
                    [yy, xdir] = log_trans_(y0,xdir0);
                    if isempty(dir(xdir))
                        mkdir(xdir)
                    end
                    silog(:,js) = redform_private(x0, yy, pshape, pd, iload, pnames, namendo, namexo, xdir, options_gsa_);
                end
                
                if isempty(threshold) && ~options_.nograph,
                    figure(hfig)
                    subplot(3,3,iplo),
                    if ilog,
                        [saso, iso] = sort(-silog(:,js));
                        bar([silog(iso(1:min(np,10)),js)])
                        logflag='log';
                    else
                        [saso, iso] = sort(-si(:,js));
                        bar(si(iso(1:min(np,10)),js))
                        logflag='';
                    end
                    %set(gca,'xticklabel',pnames(iso(1:min(np,10)),:),'fontsize',8)
                    set(gca,'xticklabel',' ','fontsize',10)
                    set(gca,'xlim',[0.5 10.5])
                    for ip=1:min(np,10),
                        text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                    title([logflag,' ',namendo,' vs. ',namexo],'interpreter','none')
                    if iplo==9,
                        dyn_saveas(hfig,[dirname,'/',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],options_);
                        if ~options_.nodisplay
                            close(hfig);
                        end
                    end
                end
                
            end
        end
    end
    if iplo<9 && iplo>0 && ifig && ~options_.nograph,
        dyn_saveas(hfig,[dirname,'/',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],options_);
        if ~options_.nodisplay
            close(hfig);
        end
    end
    ifig=0;
    iplo=0;
    for je=1:size(anamlagendo,1)
        namlagendo=deblank(anamlagendo(je,:));
        ilagendo=strmatch(namlagendo,M_.endo_names(oo_.dr.order_var(oo_.dr.nstatic+1:oo_.dr.nstatic+nsok),:),'exact');
        disp(' ')
        disp(['[', namendo,' vs. lagged ',namlagendo,']'])
        
        if ~isempty(ilagendo),
            %y0=squeeze(T(iendo,ilagendo,istable));
            y0=squeeze(T(iendo,ilagendo,:));
            if (max(y0)-min(y0))>1.e-10,
                if mod(iplo,9)==0 && isempty(threshold) && ~options_.nograph,
                    ifig=ifig+1;
                    hfig = dyn_figure(options_,'name',[namendo,' vs. lags ',int2str(ifig)]);
                    iplo=0;
                end
                iplo=iplo+1;
                js=js+1;
                xdir0 = [adir,'/',namendo,'_vs_', namlagendo];
                if ilog==0,
                    if isempty(threshold)
                        if isempty(dir(xdir0))
                            mkdir(xdir0)
                        end
                        si(:,js) = redform_private(x0, y0, pshape, pd, iload, pnames, namendo, namlagendo, xdir0, options_gsa_);
                    else
                        iy=find( (y0>threshold(1)) & (y0<threshold(2)));
                        iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
                        xdir = [xdir0,'_cut'];
                        if isempty(dir(xdir))
                            mkdir(xdir)
                        end
                        %           if ~isempty(iy)
                        %           si(:,js) = redform_private(x0(iy,:), y0(iy), pshape, pd, iload, pnames, namendo, namlagendo, xdir, options_gsa_);
                        %           end
                        if ~options_.nograph,
                            hf=dyn_figure(options_); hist(y0,30), title([namendo,' vs. ', namlagendo])
                            dyn_saveas(hf,[xdir,'/', namendo,'_vs_', namlagendo],options_);
                            if ~options_.nodisplay
                                close(hf);
                            end
                        end
                        if ~isempty(iy) && ~isempty(iyc),
                            delete([xdir, '/*cut*.*'])
                            [proba, dproba] = stab_map_1(x0, iy, iyc, 'cut',0);
                            %           indsmirnov = find(dproba>ksstat);
                            indsmirnov = find(proba<pvalue_ks);
                            for jp=1:length(indsmirnov),
                                disp([M_.param_names(estim_params_.param_vals(indsmirnov(jp),1),:),'   d-stat = ', num2str(dproba(indsmirnov(jp)),'%1.3f'),'   p-value = ', num2str(proba(indsmirnov(jp)),'%1.3f')])
                            end
                            disp(' ');
                            stab_map_1(x0, iy, iyc, 'cut',1,indsmirnov,xdir);
                            stab_map_2(x0(iy,:),alpha2,pvalue_corr,'cut',xdir)
                            stab_map_2(x0(iyc,:),alpha2,pvalue_corr,'trim',xdir)
                            lpmat=x0(iy,:);
                            lpmat0=xx0(iy,:);
                            istable=[1:length(iy)];
                            save([xdir,filesep,'threshold.mat'],'lpmat','lpmat0','istable')
                            clear lpmat lpmat0 istable
                            
                        end
                    end
                else
                    [yy, xdir] = log_trans_(y0,xdir0);
                    if isempty(dir(xdir))
                        mkdir(xdir)
                    end
                    silog(:,js) = redform_private(x0, yy, pshape, pd, iload, pnames, namendo, namlagendo, xdir, options_gsa_);
                end
                
                if isempty(threshold) && ~options_.nograph
                    figure(hfig),
                    subplot(3,3,iplo),
                    if ilog,
                        [saso, iso] = sort(-silog(:,js));
                        bar([silog(iso(1:min(np,10)),js)])
                        logflag='log';
                    else
                        [saso, iso] = sort(-si(:,js));
                        bar(si(iso(1:min(np,10)),js))
                        logflag='';
                    end
                    %set(gca,'xticklabel',pnames(iso(1:min(np,10)),:),'fontsize',8)
                    set(gca,'xticklabel',' ','fontsize',10)
                    set(gca,'xlim',[0.5 10.5])
                    for ip=1:min(np,10),
                        text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                    title([logflag,' ',namendo,' vs. ',namlagendo,'(-1)'],'interpreter','none')
                    if iplo==9,
                        dyn_saveas(hfig,[dirname,'/',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)],options_);
                        if ~options_.nodisplay
                            close(hfig);
                        end
                    end
                end
                
            end
        end
    end
    if iplo<9 && iplo>0 && ifig && ~options_.nograph,
        dyn_saveas(hfig,[dirname,'/',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)],options_);
        if ~options_.nodisplay
            close(hfig);
        end
    end
end

if isempty(threshold) && ~options_.nograph,
    if ilog==0,
        hfig=dyn_figure(options_); %bar(si)
        % boxplot(si','whis',10,'symbol','r.')
        myboxplot(si',[],'.',[],10)
        xlabel(' ')
        set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
        set(gca,'xlim',[0.5 np+0.5])
        set(gca,'ylim',[0 1])
        set(gca,'position',[0.13 0.2 0.775 0.7])
        for ip=1:np,
            text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title('Reduced form GSA')
        dyn_saveas(hfig,[dirname,'/',M_.fname,'_redform_gsa'],options_);
        
    else
        hfig=dyn_figure(options_); %bar(silog)
        % boxplot(silog','whis',10,'symbol','r.')
        myboxplot(silog',[],'.',[],10)
        set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
        xlabel(' ')
        set(gca,'xlim',[0.5 np+0.5])
        set(gca,'ylim',[0 1])
        set(gca,'position',[0.13 0.2 0.775 0.7])
        for ip=1:np,
            text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title('Reduced form GSA - Log-transformed elements')
        dyn_saveas(hfig,[dirname,'/',M_.fname,'_redform_gsa_log'],options_);
        
    end
end

function si  = redform_private(x0, y0, pshape, pd, iload, pnames, namy, namx, xdir, opt_gsa)
global bayestopt_ options_

% opt_gsa=options_.opt_gsa;
np=size(x0,2);
x00=x0;
if opt_gsa.prior_range,
    for j=1:np,
        x0(:,j)=(x0(:,j)-pd(j,3))./(pd(j,4)-pd(j,3));
    end
else
    x0=priorcdf(x0,pshape, pd(:,1), pd(:,2), pd(:,3), pd(:,4));
end

fname=[xdir,'/map'];
if iload==0,
    if isempty(dir(xdir))
        mkdir(xdir)
    end
    if ~options_.nograph,
        hfig=dyn_figure(options_); hist(y0,30), title([namy,' vs. ', namx])
        dyn_saveas(hfig,[xdir,'/', namy,'_vs_', namx],options_);
        if ~options_.nodisplay
            close(hfig);
        end
    end
    %   gsa_ = gsa_sdp_dyn(y0, x0, -2, [],[],[],1,fname, pnames);
    nrun=length(y0);
    nest=min(250,nrun);
    nfit=min(1000,nrun);
    %   dotheplots = (nfit<=nest);
    gsa_ = gsa_sdp(y0(1:nest), x0(1:nest,:), 2, [],[-1 -1 -1 -1 -1 0],[],0,[fname,'_est'], pnames);
    if nfit>nest,
        gsa_ = gsa_sdp(y0(1:nfit), x0(1:nfit,:), -2, gsa_.nvr*nest^3/nfit^3,[-1 -1 -1 -1 -1 0],[],0,fname, pnames);
    end
    save([fname,'.mat'],'gsa_')
    [sidum, iii]=sort(-gsa_.si);
    gsa_.x0=x00(1:nfit,:);
    if ~options_.nograph,
        hfig=gsa_sdp_plot(gsa_,fname,pnames,iii(1:min(12,np)));
        close(hfig);
    end
    gsa_.x0=x0(1:nfit,:);
    %   copyfile([fname,'_est.mat'],[fname,'.mat'])
    if ~options_.nograph,
        hfig=dyn_figure(options_);
        plot(y0(1:nfit),[gsa_.fit y0(1:nfit)],'.'),
        title([namy,' vs. ', namx,' fit'])
        dyn_saveas(hfig,[xdir,'/', namy,'_vs_', namx,'_fit'],options_);
        if ~options_.nodisplay
            close(hfig);
        end
        if nfit<nrun,
            npred=[nfit+1:nrun];
            yf = ss_anova_fcast(x0(npred,:), gsa_);
            hfig=dyn_figure(options_);
            plot(y0(npred),[yf y0(npred)],'.'),
            title([namy,' vs. ', namx,' pred'])
            dyn_saveas(hfig,[xdir,'/', namy,'_vs_', namx,'_pred'],options_);
            if ~options_.nodisplay
                close(hfig);
            end
        end
        
    end
else
    %   gsa_ = gsa_sdp_dyn(y0, x0, 0, [],[],[],0,fname, pnames);
    gsa_ = gsa_sdp(y0, x0, 0, [],[],[],0,fname, pnames);
    if ~options_.nograph,
        yf = ss_anova_fcast(x0, gsa_);
        hfig=dyn_figure(options_);
        plot(y0,[yf y0],'.'),
        title([namy,' vs. ', namx,' pred'])
        dyn_saveas(hfig,[xdir,'/', namy,'_vs_', namx,'_pred'],options_);
        if ~options_.nodisplay
            close(hfig);
        end
    end
end
% si = gsa_.multivariate.si;
si = gsa_.si;

