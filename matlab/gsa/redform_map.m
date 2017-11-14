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
% marco.ratto@ec.europa.eu

% Copyright (C) 2012-2016 European Commission
% Copyright (C) 2012-2017 Dynare Team
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

np = estim_params_.np;
nshock = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn;
pnames=cell(np,1);
pnames_tex=cell(np,1);
for jj=1:np
    if options_.TeX
        [param_name_temp, param_name_tex_temp]= get_the_name(nshock+jj,options_.TeX,M_,estim_params_,options_);
        pnames_tex{jj,1} = strrep(param_name_tex_temp,'$','');
        pnames{jj,1} = param_name_temp;
    else
        param_name_temp = get_the_name(nshock+jj,options_.TeX,M_,estim_params_,options_);
        pnames{jj,1} = param_name_temp;
    end
end

fname_ = M_.fname;

bounds = prior_bounds(bayestopt_, options_.prior_trunc);

if nargin==0
    dirname='';
end

if pprior
    load([dirname,filesep,M_.fname,'_prior'],'lpmat', 'lpmat0', 'istable','T');
    adir=[dirname filesep 'redform_prior'];
    type = 'prior';
else
    load([dirname,filesep,M_.fname,'_mc'],'lpmat', 'lpmat0', 'istable','T');
    adir=[dirname filesep 'redform_mc'];
    type = 'mc';
end
options_mcf.pvalue_ks = options_gsa_.ksstat_redform;
options_mcf.pvalue_corr = options_gsa_.alpha2_redform;
options_mcf.alpha2 = options_gsa_.alpha2_redform;
options_mcf.param_names = char(pnames);
if options_.TeX
    options_mcf.param_names_tex=char(pnames_tex);
end
options_mcf.fname_ = M_.fname;
options_mcf.OutputDirectoryName = adir;

if ~exist('T')
    stab_map_(dirname,options_gsa_);
    if pprior
        load([dirname,filesep,M_.fname,'_prior'],'T');
    else
        load([dirname,filesep,M_.fname,'_mc'],'T');
    end
    if ~exist('T')
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
if isempty(lpmat0)
    xx0=[];
    nshocks=0;
else
    xx0=lpmat0(istable,:);
    nshocks=size(xx0,2);
end
[kn, np]=size(x0);
offset = length(bayestopt_.pshape)-np;
if options_gsa_.prior_range
    pshape=5*(ones(np,1));
    pd =  [NaN(np,1) NaN(np,1) bounds.lb(offset+1:end) bounds.ub(offset+1:end)];
else
    pshape = bayestopt_.pshape(offset+1:end);
    pd =  [bayestopt_.p6(offset+1:end) bayestopt_.p7(offset+1:end) bayestopt_.p3(offset+1:end) bayestopt_.p4(offset+1:end)];
end

options_map.param_names = pnames;
if options_.TeX
    options_map.param_names_tex = pnames_tex;
end
options_map.fname_ = M_.fname;
options_map.OutputDirectoryName = adir;
options_map.iload = iload;
options_map.log_trans = ilog;
options_map.prior_range = options_gsa_.prior_range;
options_map.pshape = pshape;
options_map.pd = pd;

nsok = length(find(M_.lead_lag_incidence(M_.maximum_lag,:)));
lpmat=[];
lpmat0=[];
js=0;
for j=1:size(anamendo,1)
    namendo=deblank(anamendo(j,:));
    iendo=strmatch(namendo,M_.endo_names(oo_.dr.order_var,:),'exact');
    ifig=0;
    iplo=0;
    for jx=1:size(anamexo,1)
        namexo=deblank(anamexo(jx,:));
        iexo=strmatch(namexo,M_.exo_names,'exact');
        skipline()
        disp(['[', namendo,' vs ',namexo,']'])


        if ~isempty(iexo)
            %y0=squeeze(T(iendo,iexo+nspred,istable));
            y0=squeeze(T(iendo,iexo+nspred,:));
            if (max(y0)-min(y0))>1.e-10
                if mod(iplo,9)==0 && isempty(threshold) && ~options_.nograph
                    ifig=ifig+1;
                    hfig = dyn_figure(options_.nodisplay,'name',['Reduced Form Mapping: ', namendo,' vs shocks ',int2str(ifig)]);
                    iplo=0;
                end
                iplo=iplo+1;
                js=js+1;
                xdir0 = [adir,filesep,namendo,'_vs_', namexo];
                if ilog==0 || ~isempty(threshold)
                    if isempty(threshold)
                        if isempty(dir(xdir0))
                            mkdir(xdir0)
                        end
                        atitle0=['Reduced Form Mapping (ANOVA) for ',namendo,' vs ', namexo];
                        aname=[type '_' namendo '_vs_' namexo];
                        atitle=[type ' Reduced Form Mapping (ANOVA): Parameter(s) driving ',namendo,' vs ',namexo];
                        options_map.amap_name = aname;
                        options_map.amap_title = atitle;
                        options_map.figtitle = atitle0;
                        options_map.title = [namendo,' vs ', namexo];
                        options_map.OutputDirectoryName = xdir0;
                        si(:,js) = redform_private(x0, y0, options_map, options_);
                    else
                        iy=find( (y0>threshold(1)) & (y0<threshold(2)));
                        iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
                        xdir = [xdir0,'_threshold'];
                        if isempty(dir(xdir))
                            mkdir(xdir)
                        end
                        if ~options_.nograph
                            hf=dyn_figure(options_.nodisplay,'name',['Reduced Form Mapping (Monte Carlo Filtering): ',namendo,' vs ', namexo]);
                            hc = cumplot(y0);
                            a=axis; delete(hc);
                            %     hist(mat_moment{ij}),
                            x1val=max(threshold(1),a(1));
                            x2val=min(threshold(2),a(2));
                            hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'b');
                            set(hp,'FaceColor', [0.7 0.8 1])
                            hold all,
                            hc = cumplot(y0);
                            set(hc,'color','k','linewidth',2)
                            hold off,
                            title([namendo,' vs ', namexo ' - threshold [' num2str(threshold(1)) ' ' num2str(threshold(2)) ']'],'interpreter','none')
                            dyn_saveas(hf,[xdir,filesep, fname_ '_' type '_' namendo,'_vs_', namexo],options_.nodisplay,options_.graph_format);
                            create_TeX_loader(options_,[xdir,filesep, fname_ '_' type '_' namendo,'_vs_', namexo],['Reduced Form Mapping (Monte Carlo Filtering): ',strrep(namendo,'_','\_'),' vs ', strrep(namexo,'_','\_')],[type '_' namendo,'_vs_', namexo])
                        end
                        si(:,js) = NaN(np,1);
                        delete([xdir, '/*threshold*.*'])

                        atitle0=['Reduced Form Mapping (Monte Carlo Filtering) for ',namendo,' vs ', namexo];
                        aname=[type '_' namendo '_vs_' namexo '_threshold'];
                        atitle=[type ' Reduced Form Mapping (Monte Carlo Filtering): Parameter(s) driving ',namendo,' vs ',namexo];
                        options_mcf.amcf_name = aname;
                        options_mcf.amcf_title = atitle;
                        options_mcf.beha_title = 'inside threshold';
                        options_mcf.nobeha_title = 'outside threshold';
                        options_mcf.title = atitle0;
                        options_mcf.OutputDirectoryName = xdir;
                        if ~isempty(iy) && ~isempty(iyc)
                            fprintf(['%4.1f%% of the ',type,' support matches ',atitle0,'\n'],length(iy)/length(y0)*100)
                            icheck = mcf_analysis(x0, iy, iyc, options_mcf, options_);

                            lpmat=x0(iy,:);
                            if nshocks
                                lpmat0=xx0(iy,:);
                            end
                            istable=[1:length(iy)];
                            save([xdir,filesep, fname_ '_' type '_' namendo,'_vs_', namexo '_threshold' ],'lpmat','lpmat0','istable','y0','x0','xx0','iy','iyc')
                            lpmat=[]; lpmat0=[]; istable=[];
                            if length(iy)<=10 || length(iyc)<=10
                                icheck = [];  % do the generic plot in any case
                            end
                        else
                            icheck=[];
                        end
                        if isempty(icheck)
                            if length(iy)<=10 
                                if isempty(iy)
                                    disp(['There are NO MC samples in the desired range [' num2str(threshold) ']!'])
                                else
                                    disp(['There are TOO FEW (<=10) MC samples  in the desired range [' num2str(threshold) ']!'])
                                end
                            elseif length(iyc)<=10
                                if isempty(iyc)
                                    disp(['ALL MC samples are in the desired range [' num2str(threshold) ']!'])
                                else
                                    disp(['Almost ALL MC samples are in the desired range [' num2str(threshold) ']!'])
                                    disp('There are TOO FEW (<=10) MC samples OUTSIDE the desired range!')
                                end
                            end
                            atitle0=['Monte Carlo Filtering for ',namendo,' vs ', namexo];
                            options_mcf.title = atitle0;
                            indmcf = redform_mcf(y0, x0, options_mcf, options_);

                        end
                    end
                else
                    [yy, xdir] = log_trans_(y0,xdir0);
                    atitle0=['Reduced Form Mapping (ANOVA) for log-transformed ',namendo,' vs ', namexo];
                    aname=[type '_' namendo '_vs_' namexo];
                    atitle=[type ' Reduced Form Mapping (ANOVA): Parameter(s) driving ',namendo,' vs ',namexo];
                    options_map.amap_name = aname;
                    options_map.amap_title = atitle;
                    options_map.figtitle = atitle0;
                    options_map.title = ['log(' namendo ' vs ' namexo ')'];
                    options_map.OutputDirectoryName = xdir0;
                    silog(:,js) = redform_private(x0, y0, options_map, options_);
                end

                if isempty(threshold) && ~options_.nograph
                    figure(hfig)
                    subplot(3,3,iplo),
                    if ilog
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
                    for ip=1:min(np,10)
                        text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                    title([logflag,' ',namendo,' vs ',namexo],'interpreter','none')
                    if iplo==9
                        dyn_saveas(hfig,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],options_.nodisplay,options_.graph_format);
                        create_TeX_loader(options_,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],[logflag,' ',strrep(namendo,'_','\_'),' vs ',strrep(namexo,'_','\_')],['redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],1)
                    end
                end

            else
               disp(['This entry in the shock matrix is CONSTANT = ' num2str(mean(y0),3)])
            end
        end
    end
    if iplo<9 && iplo>0 && ifig && ~options_.nograph
        dyn_saveas(hfig,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],[logflag,' ',strrep(namendo,'_','\_'),' vs ',strrep(namexo,'_','\_')],['redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)],options_.figures.textwidth*min(iplo/3,1))
    end
    ifig=0;
    iplo=0;
    for je=1:size(anamlagendo,1)
        namlagendo=deblank(anamlagendo(je,:));
        ilagendo=strmatch(namlagendo,M_.endo_names(oo_.dr.order_var(M_.nstatic+1:M_.nstatic+nsok),:),'exact');
        skipline()
        disp(['[', namendo,' vs lagged ',namlagendo,']'])

        if ~isempty(ilagendo)
            %y0=squeeze(T(iendo,ilagendo,istable));
            y0=squeeze(T(iendo,ilagendo,:));
            if (max(y0)-min(y0))>1.e-10
                if mod(iplo,9)==0 && isempty(threshold) && ~options_.nograph
                    ifig=ifig+1;
                    hfig = dyn_figure(options_.nodisplay,'name',['Reduced Form Mapping: ' namendo,' vs lags ',int2str(ifig)]);
                    iplo=0;
                end
                iplo=iplo+1;
                js=js+1;
                xdir0 = [adir,filesep,namendo,'_vs_', namlagendo];
                if ilog==0 || ~isempty(threshold)
                    if isempty(threshold)
                        if isempty(dir(xdir0))
                            mkdir(xdir0)
                        end
                        atitle0=['Reduced Form Mapping (ANOVA) for ',namendo,' vs ', namlagendo];
                        aname=[type '_' namendo '_vs_' namlagendo];
                        atitle=[type ' Reduced Form Mapping (ANOVA): Parameter(s) driving ',namendo,' vs ',namlagendo];
                        options_map.amap_name = aname;
                        options_map.amap_title = atitle;
                        options_map.figtitle = atitle0;
                        options_map.title = [namendo,' vs ', namlagendo];
                        options_map.OutputDirectoryName = xdir0;
                        si(:,js) = redform_private(x0, y0, options_map, options_);
                    else
                        iy=find( (y0>threshold(1)) & (y0<threshold(2)));
                        iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
                        xdir = [xdir0,'_threshold'];
                        if isempty(dir(xdir))
                            mkdir(xdir)
                        end
                        if ~options_.nograph
                            hf=dyn_figure(options_.nodisplay,'name',['Reduced Form Mapping (Monte Carlo Filtering): ',namendo,' vs lagged ', namlagendo]);
                            hc = cumplot(y0);
                            a=axis; delete(hc);
                            %     hist(mat_moment{ij}),
                            x1val=max(threshold(1),a(1));
                            x2val=min(threshold(2),a(2));
                            hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'b');
                            set(hp,'FaceColor', [0.7 0.8 1])
                            hold all,
                            hc = cumplot(y0);
                            set(hc,'color','k','linewidth',2)
                            hold off,
                            title([namendo,' vs lagged ', namlagendo ' - threshold [' num2str(threshold(1)) ' ' num2str(threshold(2)) ']'],'interpreter','none')
                            dyn_saveas(hf,[xdir,filesep, fname_ '_' type '_' namendo,'_vs_', namlagendo],options_.nodisplay,options_.graph_format);
                            create_TeX_loader(options_,[xdir,filesep, fname_ '_' type '_' namendo,'_vs_', namlagendo],['Reduced Form Mapping (Monte Carlo Filtering): ',strrep(namendo,'_','\_'),' vs lagged ', strrep(namlagendo,'_','\_')],[type '_' namendo,'_vs_', namlagendo],1)
                        end

                        delete([xdir, '/*threshold*.*'])

                        atitle0=['Reduced Form Mapping (Monte Carlo Filtering) for ',namendo,' vs ', namlagendo];
                        aname=[type '_' namendo '_vs_' namlagendo '_threshold'];
                        atitle=[type ' Reduced Form Mapping (Monte Carlo Filtering): Parameter(s) driving ',namendo,' vs ',namlagendo];
                        options_mcf.amcf_name = aname;
                        options_mcf.amcf_title = atitle;
                        options_mcf.beha_title = 'inside threshold';
                        options_mcf.nobeha_title = 'outside threshold';
                        options_mcf.title = atitle0;
                        options_mcf.OutputDirectoryName = xdir;
                        if ~isempty(iy) && ~isempty(iyc)

                            fprintf(['%4.1f%% of the ',type,' support matches ',atitle0,'\n'],length(iy)/length(y0)*100)
                            icheck = mcf_analysis(x0, iy, iyc, options_mcf, options_);

                            lpmat=x0(iy,:);
                            if nshocks
                                lpmat0=xx0(iy,:);
                            end
                            istable=[1:length(iy)];
                            save([xdir,filesep, fname_ '_' type '_' namendo,'_vs_', namlagendo '_threshold' ],'lpmat','lpmat0','istable','y0','x0','xx0','iy','iyc')
                            lpmat=[]; lpmat0=[]; istable=[];
                            if length(iy)<=10 || length(iyc)<=10,
                                icheck = [];  % do the generic plot in any case
                            end

                        else
                            icheck = [];
                        end
                        if isempty(icheck)
                            if length(iy)<=10 
                                if isempty(iy)
                                    disp(['There are NO MC samples in the desired range [' num2str(threshold) ']!'])
                                else
                                    disp(['There are TOO FEW (<=10) MC samples  in the desired range [' num2str(threshold) ']!'])
                                end
                            elseif length(iyc)<=10
                                if isempty(iyc)
                                    disp(['ALL MC samples are in the desired range [' num2str(threshold) ']!'])
                                else
                                    disp(['Almost ALL MC samples are in the desired range [' num2str(threshold) ']!'])
                                    disp('There are TOO FEW (<=10) MC samples OUTSIDE the desired range!')
                                end
                            end
                            atitle0=['Monte Carlo Filtering for ',namendo,' vs ', namlagendo];
                            options_mcf.title = atitle0;
                            indmcf = redform_mcf(y0, x0, options_mcf, options_);
                        end
                    end
                else
                    [yy, xdir] = log_trans_(y0,xdir0);
                    atitle0=['Reduced Form Mapping (ANOVA) for log-transformed ',namendo,' vs ', namlagendo];
                    aname=[type '_' namendo '_vs_' namlagendo];
                    atitle=[type ' Reduced Form Mapping (ANOVA): Parameter(s) driving ',namendo,' vs ',namlagendo];
                    options_map.amap_name = aname;
                    options_map.amap_title = atitle;
                    options_map.figtitle = atitle0;
                    options_map.title = ['log(' namendo ' vs ' namlagendo ')'];
                    options_map.OutputDirectoryName = xdir0;
                    silog(:,js) = redform_private(x0, y0, options_map, options_);
                end

                if isempty(threshold) && ~options_.nograph
                    figure(hfig),
                    subplot(3,3,iplo),
                    if ilog
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
                    for ip=1:min(np,10)
                        text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
                    end
                    title([logflag,' ',namendo,' vs ',namlagendo,'(-1)'],'interpreter','none')
                    if iplo==9
                        dyn_saveas(hfig,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)],options_.nodisplay,options_.graph_format);
                        create_TeX_loader(options_,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)],[logflag,' ',strrep(namendo,'_','\_'),' vs ',strrep(namlagendo,'_','\_'),'(-1)'],['redform_', namendo,'_vs_lags_',logflag,':',num2str(ifig)],1)
                    end
                end

            else
               disp(['This entry in the transition matrix is CONSTANT = ' num2str(mean(y0),3)])
            end
        end
    end
    if iplo<9 && iplo>0 && ifig && ~options_.nograph
        dyn_saveas(hfig,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[dirname,filesep,M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)],[logflag,' ',strrep(namendo,'_','\_'),' vs ',strrep(namlagendo,'_','\_'),'(-1)'],['redform_', namendo,'_vs_lags_',logflag,':',num2str(ifig)],options_.figures.textwidth*min(iplo/3,1));
    end
end

if isempty(threshold) && ~options_.nograph
    if ilog==0
        hfig=dyn_figure(options_.nodisplay,'name','Reduced Form GSA'); %bar(si)
                                                                       % boxplot(si','whis',10,'symbol','r.')
        myboxplot(si',[],'.',[],10)
        xlabel(' ')
        set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
        set(gca,'xlim',[0.5 np+0.5])
        set(gca,'ylim',[0 1])
        set(gca,'position',[0.13 0.2 0.775 0.7])
        for ip=1:np
            text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title('Reduced form GSA')
        dyn_saveas(hfig,[dirname,filesep,M_.fname,'_redform_gsa'],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[dirname,filesep,M_.fname,'_redform_gsa'],'Reduced Form GSA','redform_gsa')

    else
        hfig=dyn_figure(options_.nodisplay,'name','Reduced Form GSA'); %bar(silog)
                                                                       % boxplot(silog','whis',10,'symbol','r.')
        myboxplot(silog',[],'.',[],10)
        set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
        xlabel(' ')
        set(gca,'xlim',[0.5 np+0.5])
        set(gca,'ylim',[0 1])
        set(gca,'position',[0.13 0.2 0.775 0.7])
        for ip=1:np
            text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title('Reduced form GSA - Log-transformed elements')
        dyn_saveas(hfig,[dirname,filesep,M_.fname,'_redform_gsa_log'],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[dirname,filesep,M_.fname,'_redform_gsa_log'],'Reduced form GSA - Log-transformed elements','redform_gsa_log')

    end
end

function si  = redform_private(x0, y0, options_map, options_)

np=size(x0,2);
x00=x0;
ilog = options_map.log_trans;
iload = options_map.iload;
pnames = options_map.param_names;
pd = options_map.pd;
pshape = options_map.pshape;
xdir = options_map.OutputDirectoryName;
if options_map.prior_range
    for j=1:np
        x0(:,j)=(x0(:,j)-pd(j,3))./(pd(j,4)-pd(j,3));
    end
else
    x0=priorcdf(x0,pshape, pd(:,1), pd(:,2), pd(:,3), pd(:,4));
end

if ilog
    fname=[xdir filesep options_map.fname_ '_' options_map.amap_name '_log'];
else
    fname=[xdir filesep options_map.fname_ '_' options_map.amap_name];
end
if iload==0
    if isempty(dir(xdir))
        mkdir(xdir)
    end
    nrun=length(y0);
    nest=max(50,nrun/2);
    nest=min(250,nest);
    nfit=min(1000,nrun);
    %   dotheplots = (nfit<=nest);
    %     gsa_ = gsa_sdp(y0(1:nest), x0(1:nest,:), 2, [],[-1 -1 -1 -1 -1 0],[],0,[fname,'_est'], pnames);
    [ys,is] = sort(y0);
    istep = ceil(nrun/nest);
    if istep>1
        iest = is(floor(istep/2):istep:end);
        nest = length(iest);
        irest = is(setdiff([1:nrun],[floor(istep/2):istep:nrun]));
        istep = ceil(length(irest)/(nfit-nest));
        ifit = union(iest, irest(1:istep:end));
    else
        warning('the number of samples is too small for ANOVA estimation')
        si=nan(np,1);
        return
    end
    if ~ismember(irest(end),ifit)
        ifit = union(ifit, irest(end));
    end
    nfit=length(ifit);
    %     ifit = union(iest, irest(randperm(nrun-nest,nfit-nest)));
    %     ifit = iest;
    %     nfit=nest;
    ipred = setdiff([1:nrun],ifit);

    if ilog
        [y1, tmp, isig, lam] = log_trans_(y0(iest));
        y1 = log(y0*isig+lam);
    end
    if ~options_.nograph
        hfig=dyn_figure(options_.nodisplay,'name',options_map.figtitle);
        subplot(221)
        if ilog
            hist(y1,30)
        else
            hist(y0,30)
        end
        title(options_map.title,'interpreter','none')
        subplot(222)
        if ilog
            hc = cumplot(y1);
        else
            hc = cumplot(y0);
        end
        set(hc,'color','k','linewidth',2)
        title([options_map.title ' CDF'],'interpreter','none')
    end

    gsa0 = ss_anova(y0(iest), x0(iest,:), 1);
    if ilog
        [gsa22, gsa1, gsax] = ss_anova_log(y1(iest), x0(iest,:), isig, lam, gsa0);
    end
    %     if (gsa1.out.bic-gsa0.out.bic) < 10,
    %         y00=y0;
    %         gsa00=gsa0;
    %         gsa0=gsa1;
    %         y0=y1;
    %         ilog=1;
    %     end
    if nfit>nest
        %         gsa_ = gsa_sdp(y0(1:nfit), x0(1:nfit,:), -2, gsa_.nvr*nest^3/nfit^3,[-1 -1 -1 -1 -1 0],[],0,fname, pnames);
        nvr =  gsa0.nvr*nest^3/nfit^3;
        nvr(gsa0.stat<2) = gsa0.nvr(gsa0.stat<2)*nest^5/nfit^5;
        gsa_ = ss_anova(y0(ifit), x0(ifit,:), 1, 0, 2, nvr);
        if ilog
            gsa0 = gsa_;
            nvr1 =  gsa1.nvr*nest^3/nfit^3;
            nvr1(gsa1.stat<2) = gsa1.nvr(gsa1.stat<2)*nest^5/nfit^5;
            nvrx =  gsax.nvr*nest^3/nfit^3;
            nvrx(gsax.stat<2) = gsax.nvr(gsax.stat<2)*nest^5/nfit^5;
            [gsa22, gsa1, gsax] = ss_anova_log(y1(ifit), x0(ifit,:), isig, lam, gsa0, [nvr1' nvrx']);
            %         gsa1 = ss_anova(y1(ifit), x0(ifit,:), 1, 0, 2, nvr);
            %         gsa2=gsa1;
            %         gsa2.y = gsa0.y;
            %         gsa2.fit = (exp(gsa1.fit)-lam)*isig;
            %         gsa2.f0 = mean(gsa2.fit);
            %         gsa2.out.SSE = sum((gsa2.fit-gsa2.y).^2);
            %         gsa2.out.bic = gsa2.out.bic-nest*log(gsa1.out.SSE)+nest*log(gsa2.out.SSE);
            %         gsa2.r2 = 1-cov(gsa2.fit-gsa2.y)/cov(gsa2.y);
            %         for j=1:np,
            %             gsa2.fs(:,j) = exp(gsa1.fs(:,j)).*mean(exp(gsa1.fit-gsa1.f(:,j)))*isig-lam*isig-gsa2.f0;
            %             gsa2.f(:,j) = exp(gsa1.f(:,j)).*mean(exp(gsa1.fit-gsa1.f(:,j)))*isig-lam*isig-gsa2.f0;
            %             gsa2.si(j) = var(gsa2.f(:,j))/var(gsa2.y);
            %         end
            %         nvr =  gsax.nvr*nest^3/nfit^3;
            %         nvr(gsax.stat<2) = gsax.nvr(gsax.stat<2)*nest^5/nfit^5;
            %         gsax = ss_anova([gsa2.y-gsa2.fit], x0(ifit,:), 1, 0, 2, nvr);
            %         gsa22=gsa2;
            %         gsa22.fit = gsa2.fit+gsax.fit;
            %         gsa22.f0 = mean(gsa22.fit);
            %         gsa22.out.SSE = sum((gsa22.fit-gsa22.y).^2);
            %         gsa22.out.bic = nest*log(gsa22.out.SSE/nest) + (gsax.out.df+gsa2.out.df-1)*log(nest);
            %         gsa22.r2 = 1-sum((gsa22.fit-gsa22.y).^2)/sum((gsa22.y-mean(gsa22.y)).^2);
            %         for j=1:np,
            %             gsa22.fs(:,j) = gsa2.fs(:,j)+gsax.fs(:,j);
            %             gsa22.f(:,j) = gsa2.f(:,j)+gsax.f(:,j);
            %             gsa22.si(j) = var(gsa22.f(:,j))/var(gsa22.y);
            %         end
            gsa_ = gsa22;
        end
    else
        if ilog
            gsa_ = gsa22;
        else
            gsa_ = gsa0;
        end
    end
    save([fname,'_map.mat'],'gsa_')
    [sidum, iii]=sort(-gsa_.si);
    gsa_.x0=x00(ifit,:);
    if ~options_.nograph
        hmap=gsa_sdp_plot(gsa_,[fname '_map'],pnames,iii(1:min(12,np)));
        set(hmap,'name',options_map.amap_title);
    end
    gsa_.x0=x0(ifit,:);
    %   copyfile([fname,'_est.mat'],[fname,'.mat'])
    if ~options_.nograph
        figure(hfig);
        subplot(223),
        plot(y0(ifit),[gsa_.fit y0(ifit)],'.'),
        r2 = gsa_.r2;
        %         if ilog,
        %             plot(y00(ifit),[log_trans_(gsa_.fit,'',isig,lam) y00(ifit)],'.'),
        %             r2 = 1 - cov(log_trans_(gsa_.fit,'',isig,lam)-y00(ifit))/cov(y00(ifit));
        %         else
        %             plot(y0(ifit),[gsa_.fit y0(ifit)],'.'),
        %             r2 = gsa_.r2;
        %         end
        title(['Learning sample fit - R2=' num2str(r2,2)],'interpreter','none')
        if nfit<nrun
            if ilog
                yf = ss_anova_fcast(x0(ipred,:), gsa1);
                yf = log_trans_(yf,'',isig,lam)+ss_anova_fcast(x0(ipred,:), gsax);
            else
                yf = ss_anova_fcast(x0(ipred,:), gsa_);
            end
            yn = y0(ipred);
            r2  = 1-cov(yf-yn)/cov(yn);
            subplot(224),
            plot(yn,[yf yn],'.'),
            title(['Out-of-sample prediction - R2=' num2str(r2,2)],'interpreter','none')
        end
        dyn_saveas(hfig,fname,options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,fname,['Out-of-sample prediction - R2=' num2str(r2,2)],'redform_gsa_log')

        if options_.nodisplay
            close(hmap);
        end
    end
else
    %   gsa_ = gsa_sdp_dyn(y0, x0, 0, [],[],[],0,fname, pnames);
    %     gsa_ = gsa_sdp(y0, x0, 0, [],[],[],0,fname, pnames);
    load([fname,'_map.mat'],'gsa_')
    if ~options_.nograph
        yf = ss_anova_fcast(x0, gsa_);
        hfig=dyn_figure(options_.nodisplay,'name',options_map.title);
        plot(y0,[yf y0],'.'),
        title([namy,' vs ', namx,' pred'],'interpreter','none')
        dyn_saveas(hfig,[fname '_pred'],options_.nodisplay,options_.graph_format);
        create_TeX_loader(options_,[fname '_pred'],options_map.title,[namy,' vs ', namx,' pred'])

    end
end
% si = gsa_.multivariate.si;
si = gsa_.si;

return

function gsa2 = log2level_map(gsa1, isig, lam)

nest=length(gsa1.y);
np = size(gsa1.x0,2);
gsa2=gsa1;
gsa2.y = log_trans_(gsa1.y,'',isig,lam);
gsa2.fit = (exp(gsa1.fit)-lam)*isig;
gsa2.f0 = mean(gsa2.fit);
gsa2.out.SSE = sum((gsa2.fit-gsa2.y).^2);
gsa2.out.bic = gsa2.out.bic-nest*log(gsa1.out.SSE)+nest*log(gsa2.out.SSE);
gsa2.r2 = 1-cov(gsa2.fit-gsa2.y)/cov(gsa2.y);
for j=1:np
    gsa2.fs(:,j) = exp(gsa1.fs(:,j)).*mean(exp(gsa1.fit-gsa1.f(:,j)))*isig-lam*isig-gsa2.f0;
    gsa2.fses(:,j) = exp(gsa1.fs(:,j)+gsa1.fses(:,j)).*mean(exp(gsa1.fit-gsa1.f(:,j)))*isig-lam*isig-gsa2.f0-gsa2.fs(:,j);
    gsa2.f(:,j) = exp(gsa1.f(:,j)).*mean(exp(gsa1.fit-gsa1.f(:,j)))*isig-lam*isig-gsa2.f0;
    gsa2.si(j) = var(gsa2.f(:,j))/var(gsa2.y);
end

return


function [gsa22, gsa1, gsax] = ss_anova_log(y,x,isig,lam,gsa0,nvrs)

[nest, np]=size(x);

if nargin==6
    gsa1 = ss_anova(y, x, 1, 0, 2, nvrs(:,1));
else
    gsa1 = ss_anova(y, x, 1);
end
gsa2 = log2level_map(gsa1, isig, lam);
if nargin >=5 && ~isempty(gsa0)
    for j=1:np
        nvr2(j) = var(diff(gsa2.fs(:,j),2));
        nvr0(j) = var(diff(gsa0.fs(:,j),2));
    end
    inda = find((gsa0.stat<2)&(gsa1.stat>2));
    inda = inda(log10(nvr0(inda)./nvr2(inda))/2<0);
    gsa1.nvr(inda)=gsa1.nvr(inda).*10.^(log10(nvr0(inda)./nvr2(inda)));
    gsa1 = ss_anova(y, x, 1, 0, 2, gsa1.nvr);
    gsa2 = log2level_map(gsa1, isig, lam);
end
if nargin==6
    gsax = ss_anova(gsa2.y-gsa2.fit, x, 1, 0, 2, nvrs(:,2));
else
    gsax = ss_anova(gsa2.y-gsa2.fit, x, 1);
end
gsa22=gsa2;
gsa22.fit = gsa2.fit+gsax.fit;
gsa22.f0 = mean(gsa22.fit);
gsa22.out.SSE = sum((gsa22.fit-gsa22.y).^2);
gsa22.out.bic = nest*log(gsa22.out.SSE/nest) + (gsax.out.df+gsa2.out.df-1)*log(nest);
gsa22.r2 = 1-sum((gsa22.fit-gsa22.y).^2)/sum((gsa22.y-mean(gsa22.y)).^2);
for j=1:np
    gsa22.fs(:,j) = gsa2.fs(:,j)+gsax.fs(:,j);
    gsa22.fses(:,j) = gsax.fses(:,j);
    gsa22.f(:,j) = gsa2.f(:,j)+gsax.f(:,j);
    gsa22.si(j) = var(gsa22.f(:,j))/var(gsa22.y);
end

return

function indmcf = redform_mcf(y0, x0, options_mcf, options_)

hfig=dyn_figure(options_.nodisplay,'name',options_mcf.amcf_title);

[post_mean, post_median, post_var, hpd_interval, post_deciles, ...
 density] = posterior_moments(y0,1,0.9);
post_deciles = [-inf; post_deciles; inf];

for jt=1:10
    indy{jt}=find( (y0>post_deciles(jt)) & (y0<=post_deciles(jt+1)));
    leg{jt}=[int2str(jt) '-dec'];
end
[proba, dproba] = stab_map_1(x0, indy{1}, indy{end}, [],0);
indmcf=find(proba<options_mcf.pvalue_ks);
if isempty(indmcf)
    [tmp,jtmp] = sort(proba,2,'ascend');
    indmcf = jtmp(1);
%     indmcf = jtmp(1:min(2,length(proba)));
end
[tmp,jtmp] = sort(proba(indmcf),2,'ascend');
indmcf = indmcf(jtmp);
nbr_par = length(indmcf);
nrow=ceil(sqrt(nbr_par+1));
ncol=nrow;
if nrow*(nrow-1)>nbr_par
    ncol=nrow-1;
end

cmap = colormap(jet(10));
for jx=1:nbr_par
    subplot(nrow,ncol,jx)
    hold off
    for jt=1:10
        h=cumplot(x0(indy{jt},indmcf(jx)));
        set(h,'color', cmap(jt,:), 'linewidth', 2)
        hold all
    end
    title(options_mcf.param_names(indmcf(jx),:),'interpreter','none')
end
hleg = legend(leg);
aa=get(hleg,'Position');
aa(1)=1-aa(3)-0.02;
aa(2)=0.02;
set(hleg,'Position',aa);
if ~isoctave
    annotation('textbox', [0.25,0.01,0.5,0.05], ...
               'String', options_mcf.title, ...
               'Color','black',...
               'FontWeight','bold',...
               'interpreter','none',...
               'horizontalalignment','center');
end

dyn_saveas(hfig,[options_mcf.OutputDirectoryName filesep options_mcf.fname_,'_',options_mcf.amcf_name],options_.nodisplay,options_.graph_format);
create_TeX_loader(options_,[options_mcf.OutputDirectoryName filesep options_mcf.fname_,'_',options_mcf.amcf_name],strrep(options_mcf.amcf_title,'_','\_'),[options_mcf.fname_,'_',options_mcf.amcf_name])

return

function []=create_TeX_loader(options_,figpath,caption,label_name,scale_factor)
if nargin<5
    scale_factor=1;
end
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([figpath '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by redform_map.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',scale_factor,strrep(figpath,'\','/'));
    fprintf(fidTeX,'\\caption{%s.}',caption);
    fprintf(fidTeX,'\\label{Fig:%s}\n',label_name);
    fprintf(fidTeX,'\\end{figure}\n\n');
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end
