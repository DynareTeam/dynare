function x0=dynare_sensitivity(options_gsa)
% Frontend to the Sensitivity Analysis Toolbox for DYNARE
%
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.

% Copyright (C) 2008-2011 Dynare Team
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

global M_ options_ oo_ bayestopt_ estim_params_

fname_ = M_.fname;
if ~isfield(M_,'dname'),
    M_.dname = M_.fname;
end
lgy_ = M_.endo_names;
x0=[];

options_gsa = set_default_option(options_gsa,'datafile',[]);
options_gsa = set_default_option(options_gsa,'rmse',0);
if isfield(options_gsa,'nograph'),
    options_.nograph=options_gsa.nograph;
end
if isfield(options_gsa,'nodisplay'),
    options_.nodisplay=options_gsa.nodisplay;
end
if isfield(options_gsa,'graph_format'),
    options_.graph_format=options_gsa.graph_format;
end
if isfield(options_gsa,'mode_file'),
    options_.mode_file=options_gsa.mode_file;
end

options_.order = 1;

if ~isempty(options_gsa.datafile) || isempty(bayestopt_) || options_gsa.rmse,
    if isempty(options_gsa.datafile) && options_gsa.rmse,
        disp('The data file and all relevant estimation options ')
        disp('[first_obs, nobs, presample, prefilter, loglinear, lik_init, kalman_algo, ...]')
        disp('must be specified for RMSE analysis!');
        error('Sensitivity anaysis error!');
    end
    
    options_.datafile = options_gsa.datafile;
    if isfield(options_gsa,'first_obs'),
        options_.first_obs=options_gsa.first_obs;
    end
    if isfield(options_gsa,'nobs'),
        options_.nobs=options_gsa.nobs;
    end
    if isfield(options_gsa,'presample'),
        options_.presample=options_gsa.presample;
    end
    if isfield(options_gsa,'prefilter'),
        options_.prefilter=options_gsa.prefilter;
    end
    if isfield(options_gsa,'loglinear'),
        options_.loglinear=options_gsa.loglinear;
    end
    if isfield(options_gsa,'lik_init'),
        options_.lik_init=options_gsa.lik_init;
    end
    if isfield(options_gsa,'kalman_algo'),
        options_.kalman_algo=options_gsa.kalman_algo;
    end
    options_.mode_compute = 0;
    options_.filtered_vars = 1;
    options_.plot_priors = 0;
%     [data,rawdata,xparam1,data_info]=dynare_estimation_init([],fname_,1);
    [dataset_,xparam1, M_, options_, oo_, estim_params_,bayestopt_]=dynare_estimation_init(M_.endo_names,fname_,1, M_, options_, oo_, estim_params_, bayestopt_);
    % computes a first linear solution to set up various variables
else
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    end    
end
[make,my,day,punk,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);

options_gsa = set_default_option(options_gsa,'identification',0);
if options_gsa.identification,
    options_gsa.redform=0;
    options_gsa = set_default_option(options_gsa,'morris',1);
    options_gsa = set_default_option(options_gsa,'trans_ident',0);
    options_gsa = set_default_option(options_gsa,'load_ident_files',0);
    options_gsa = set_default_option(options_gsa,'ar',1);
    options_gsa = set_default_option(options_gsa,'useautocorr',0);
    options_.ar = options_gsa.ar;
    if options_gsa.morris==0,
        disp('The option morris = 0 is no longer supported (Type I errors)')
        disp('This option is reset at morris = 2 (local identification analysis).')
        options_gsa.morris=2;
    end
    if options_gsa.morris==2,
        if isfield(options_,'options_ident'),
            options_.options_ident.load_ident_files = options_gsa.load_ident_files;
            options_.options_ident.useautocorr = options_gsa.useautocorr;
            options_.options_ident.ar = options_gsa.ar;            
            options_ident=options_.options_ident;
        else
            options_ident=[];
            options_ident = set_default_option(options_ident,'load_ident_files',options_gsa.load_ident_files);
            options_ident = set_default_option(options_ident,'useautocorr',options_gsa.useautocorr);
            options_ident = set_default_option(options_ident,'ar',options_gsa.ar);
            options_.options_ident = options_ident;
        end
    end
end

% map stability
options_gsa = set_default_option(options_gsa,'stab',1);
options_gsa = set_default_option(options_gsa,'redform',0);
options_gsa = set_default_option(options_gsa,'pprior',1);
options_gsa = set_default_option(options_gsa,'prior_range',1);
options_gsa = set_default_option(options_gsa,'ppost',0);
options_gsa = set_default_option(options_gsa,'neighborhood_width',0);
options_gsa = set_default_option(options_gsa,'ilptau',1);
options_gsa = set_default_option(options_gsa,'morris',0);
options_gsa = set_default_option(options_gsa,'glue',0);
options_gsa = set_default_option(options_gsa,'morris_nliv',6);
options_gsa = set_default_option(options_gsa,'morris_ntra',20);
options_gsa = set_default_option(options_gsa,'Nsam',2048);
options_gsa = set_default_option(options_gsa,'load_redform',0);
options_gsa = set_default_option(options_gsa,'load_rmse',0);
options_gsa = set_default_option(options_gsa,'load_stab',0);
options_gsa = set_default_option(options_gsa,'alpha2_stab',0);
options_gsa = set_default_option(options_gsa,'ksstat',0.1);
options_gsa = set_default_option(options_gsa,'pvalue_ks',0.001);
options_gsa = set_default_option(options_gsa,'pvalue_corr',0.001);
%options_gsa = set_default_option(options_gsa,'load_mh',0);
% REDFORM mapping
options_gsa = set_default_option(options_gsa,'logtrans_redform',0);
options_gsa = set_default_option(options_gsa,'threshold_redform',[]);
options_gsa = set_default_option(options_gsa,'ksstat_redform',0.001);
options_gsa = set_default_option(options_gsa,'alpha2_redform',0);
options_gsa = set_default_option(options_gsa,'namendo',[]);
options_gsa = set_default_option(options_gsa,'namlagendo',[]);
options_gsa = set_default_option(options_gsa,'namexo',[]);
% RMSE mapping
options_gsa = set_default_option(options_gsa,'lik_only',0);
options_gsa = set_default_option(options_gsa,'var_rmse',options_.varobs);
options_gsa = set_default_option(options_gsa,'pfilt_rmse',0.1);
options_gsa = set_default_option(options_gsa,'istart_rmse',options_.presample+1);
options_gsa = set_default_option(options_gsa,'alpha_rmse',0.001);
options_gsa = set_default_option(options_gsa,'alpha2_rmse',0);

if options_gsa.redform && options_gsa.neighborhood_width==0 && isempty(options_gsa.threshold_redform),
    options_gsa.pprior=1;
    options_gsa.ppost=0;
end

if options_gsa.morris>2,
    disp('The option morris = 3 is no longer supported')
    disp('the option is reset at morris = 1 .')
    options_gsa.morris=1;
end
   
if options_gsa.morris==1,
    if ~options_gsa.identification,
        options_gsa.redform=1;
    end
    options_gsa.pprior=1;
    options_gsa.ppost=0;
    %options_gsa.stab=1;
    options_gsa.glue=0;
    options_gsa.rmse=0;
    options_gsa.load_rmse=0;
    options_gsa.alpha2_stab=1;
    options_gsa.ksstat=1;
    options_gsa.pvalue_ks=0;
    options_gsa.pvalue_corr=0;
%     if options_gsa.morris==3,
%         options_gsa = set_default_option(options_gsa,'Nsam',256);
%         OutputDirectoryName = CheckPath('gsa/identif',M_.dname);
%     else
        OutputDirectoryName = CheckPath('gsa/screen',M_.dname);
%     end
else
    OutputDirectoryName = CheckPath('gsa',M_.dname);
end

% options_.opt_gsa = options_gsa;

if (options_gsa.load_stab || options_gsa.load_rmse || options_gsa.load_redform) ...
        && options_gsa.pprior,
    filetoload=[OutputDirectoryName '/' fname_ '_prior.mat'];
    if ~exist(filetoload),
        disp([filetoload,' not found!'])
        disp(['You asked to load a non existent analysis'])
        %options_gsa.load_stab=0;
        return,
    else
        if isempty(strmatch('bkpprior',who('-file', filetoload),'exact')),
            disp('Warning! Missing prior info for saved sample') % trap for files previous 
            disp('The saved files are generated with previous version of GSA package') % trap for files previous 
        else
            load(filetoload,'bkpprior'),
            if any(bayestopt_.pshape~=bkpprior.pshape) || ...
                    any(bayestopt_.p1~=bkpprior.p1) || ...
                    any(bayestopt_.p2~=bkpprior.p2) || ...
                    any(bayestopt_.p3(~isnan(bayestopt_.p3))~=bkpprior.p3(~isnan(bkpprior.p3))) || ...
                    any(bayestopt_.p4(~isnan(bayestopt_.p4))~=bkpprior.p4(~isnan(bkpprior.p4))),
                disp('WARNING!')
                disp('The saved sample has different priors w.r.t. to current ones!!')
                disp('')
                disp('Press ENTER to continue')
                pause;
            end
        end
    end
end

if options_gsa.stab && ~options_gsa.ppost,
    x0 = stab_map_(OutputDirectoryName,options_gsa);
    if isempty(x0),
        disp(' ')
        disp('Sensitivity computations stopped: no parameter set provided a unique solution')
        return
    end
end

% reduced form
% redform_map(namendo, namlagendo, namexo, icomp, pprior, ilog, threshold)

options_.opt_gsa = options_gsa;
if options_gsa.identification,
    map_ident_(OutputDirectoryName,options_gsa);
end

if options_gsa.redform && ~isempty(options_gsa.namendo),% ...
    %         && ~options_gsa.ppost,
    if options_gsa.ppost,
        filnam = dir([M_.dname filesep 'metropolis' filesep '*param_irf*.mat']);
        lpmat=[];
        for j=1:length(filnam),
            load ([M_.dname filesep 'metropolis' filesep M_.fname '_param_irf' int2str(j) '.mat'])
            lpmat=[lpmat; stock];
        end
        clear stock
        nshock = estim_params_.nvx;
        nshock = nshock + estim_params_.nvn;
        nshock = nshock + estim_params_.ncx;
        nshock = nshock + estim_params_.ncn;
        
        lpmat0=lpmat(:,1:nshock);
        lpmat=lpmat(:,nshock+1:end);
        istable=(1:size(lpmat,1));
        iunstable=[];
        iwrong=[];
        iindeterm=[];
        save([OutputDirectoryName filesep M_.fname '_mc.mat'],'lpmat','lpmat0','istable','iunstable','iwrong','iindeterm')
        options_gsa.load_stab=1;
        
        x0 = stab_map_(OutputDirectoryName,options_gsa);
    end
    if strmatch(':',options_gsa.namendo,'exact'),
        options_gsa.namendo=M_.endo_names(1:M_.orig_endo_nbr,:);
    end
    if strmatch(':',options_gsa.namexo,'exact'),
        options_gsa.namexo=M_.exo_names;
    end
    if strmatch(':',options_gsa.namlagendo,'exact'),
        options_gsa.namlagendo=M_.endo_names(1:M_.orig_endo_nbr,:);
    end
%     options_.opt_gsa = options_gsa;
    if options_gsa.morris==1,
        redform_screen(OutputDirectoryName,options_gsa);
    else
        % check existence of the SS_ANOVA toolbox
        if isempty(options_gsa.threshold_redform) && ...
         ~(exist('gsa_sdp','file')==6 || exist('gsa_sdp','file')==2),
            disp('Download Mapping routines at:')
            disp('http://eemc.jrc.ec.europa.eu/softwareDYNARE-Dowload.htm')
            disp(' ' )
            error('Mapping routines missing!')
        end

        redform_map(OutputDirectoryName,options_gsa);
    end
end
% RMSE mapping
% function [rmse_MC, ixx] = filt_mc_(vvarvecm, loadSA, pfilt, alpha, alpha2)
options_.opt_gsa = options_gsa;
if options_gsa.rmse,
    if ~options_gsa.ppost
        if options_gsa.pprior
            a=whos('-file',[OutputDirectoryName,'/',fname_,'_prior'],'logpo2');
        else
            a=whos('-file',[OutputDirectoryName,'/',fname_,'_mc'],'logpo2');
        end
        if exist('OCTAVE_VERSION'),
            aflag=0;
            for ja=1:length(a),
                aflag=aflag+strcmp('logpo2',a(ja).name);
            end
            if aflag==0,
                a=[];
            else
                a=1;
            end
        end
        if isempty(a),
           if options_gsa.lik_only,
               options_.smoother=0;
               options_.filter_step_ahead=[];
               options_.forecast=0;
               options_.filtered_vars=0;               
           end
%             dynare_MC([],OutputDirectoryName,data,rawdata,data_info);
        if options_gsa.pprior
            TmpDirectoryName = ([M_.dname filesep 'gsa' filesep 'prior']);
        else
            TmpDirectoryName = ([M_.dname filesep 'gsa' filesep 'mc']);
        end
        if exist(TmpDirectoryName,'dir');
            mydelete([M_.fname '_filter_step_ahead*.mat'],[TmpDirectoryName filesep]);
            mydelete([M_.fname '_inno*.mat'],[TmpDirectoryName filesep]);
            mydelete([M_.fname '_smooth*.mat'],[TmpDirectoryName filesep]);
            mydelete([M_.fname '_update*.mat'],[TmpDirectoryName filesep]);
            filparam = dir([TmpDirectoryName filesep M_.fname '_param*.mat']);
            for j=1:length(filparam),
                if isempty(strmatch([M_.fname '_param_irf'],filparam(j).name))
                    delete([TmpDirectoryName filesep filparam(j).name]);
                end
            end
            
        end
        prior_posterior_statistics('gsa',dataset_);
        if options_.bayesian_irf
            PosteriorIRF('gsa');
        end
        options_gsa.load_rmse=0;
        %   else
        %     if options_gsa.load_rmse==0,
        %       disp('You already saved a MC filter/smoother analysis ')
        %       disp('Do you want to overwrite ?')
        %       pause;
        %       if options_gsa.pprior
        %         delete([OutputDirectoryName,'/',fname_,'_prior_*.mat'])
        %       else
        %         delete([OutputDirectoryName,'/',fname_,'_mc_*.mat'])
        %       end
        %       dynare_MC([],OutputDirectoryName);
        %       options_gsa.load_rmse=0;
        %     end
            
        end
    end
    clear a;
%     filt_mc_(OutputDirectoryName,data_info);
    filt_mc_(OutputDirectoryName,options_gsa,dataset_);
end
options_.opt_gsa = options_gsa;


if options_gsa.glue,
    dr_ = oo_.dr;
    if options_gsa.ppost
        load([OutputDirectoryName,'/',fname_,'_post']);
        DirectoryName = CheckPath('metropolis',M_.dname);
    else
        if options_gsa.pprior
            load([OutputDirectoryName,'/',fname_,'_prior']);
        else
            load([OutputDirectoryName,'/',fname_,'_mc']);
        end
    end
    if ~exist('x'),
        disp(['No RMSE analysis is available for current options'])
        disp(['No GLUE file prepared'])
        return,
    end
    nruns=size(x,1);
    gend = options_.nobs;
    rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
    rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
    if options_.loglinear == 1
        rawdata = log(rawdata);
    end
    if options_.prefilter == 1
        %data = transpose(rawdata-ones(gend,1)*bayestopt_.mean_varobs);
        data = transpose(rawdata-ones(gend,1)*mean(rawdata,1));
    else
        data = transpose(rawdata);
    end
    
    Obs.data = data;
    Obs.time = [1:gend];
    Obs.num  = gend;
    for j=1:size(options_.varobs,1)
        Obs.name{j} = deblank(options_.varobs(j,:));
        vj=deblank(options_.varobs(j,:));
        
        jxj = strmatch(vj,lgy_(dr_.order_var,:),'exact');
        js = strmatch(vj,lgy_,'exact');
        if ~options_gsa.ppost
            %       y0=zeros(gend+1,nruns);
            %       nb = size(stock_filter,3);
            %       y0 = squeeze(stock_filter(:,jxj,:)) + ...
            %         kron(stock_ys(js,:),ones(size(stock_filter,1),1));
            %       Out(j).data = y0';
            %       Out(j).time = [1:size(y0,1)];
            Out(j).data = jxj;
            Out(j).time = [pwd,'/',OutputDirectoryName];
        else
            Out(j).data = jxj;
            Out(j).time = [pwd,'/',DirectoryName];
        end
        Out(j).name = vj;
        Out(j).ini  = 'yes';
        Lik(j).name = ['rmse_',vj];
        Lik(j).ini  = 'yes';
        Lik(j).isam = 1;
        Lik(j).data = rmse_MC(:,j)';
        
        if ~options_gsa.ppost
            %       y0 = squeeze(stock_smooth(:,jxj,:)) + ...
            %         kron(stock_ys(js,:),ones(size(stock_smooth,1),1));
            %       Out1(j).name = vj;
            %       Out1(j).ini  = 'yes';
            %       Out1(j).time = [1:size(y0,1)];
            %       Out1(j).data = y0';
            Out1=Out;
        else
            Out1=Out;
        end
        ismoo(j)=jxj;
        
    end
    jsmoo = size(options_.varobs,1);
    for j=1:M_.endo_nbr,
        if ~ismember(j,ismoo),
            jsmoo=jsmoo+1;
            vj=deblank(M_.endo_names(dr_.order_var(j),:));
            if ~options_gsa.ppost        
                %         y0 = squeeze(stock_smooth(:,j,:)) + ...
                %           kron(stock_ys(j,:),ones(size(stock_smooth,1),1));
                %         Out1(jsmoo).time = [1:size(y0,1)];
                %         Out1(jsmoo).data = y0';
                Out1(jsmoo).data = j;
                Out1(jsmoo).time = [pwd,'/',OutputDirectoryName];
            else
                Out1(jsmoo).data = j;
                Out1(jsmoo).time = [pwd,'/',DirectoryName];
            end
            Out1(jsmoo).name = vj;
            Out1(jsmoo).ini  = 'yes';
        end
    end
    tit(M_.exo_names_orig_ord,:) = M_.exo_names;
    for j=1:M_.exo_nbr,
        Exo(j).name = deblank(tit(j,:));    
    end
    if ~options_gsa.ppost
        Lik(size(options_.varobs,1)+1).name = 'logpo';
        Lik(size(options_.varobs,1)+1).ini  = 'yes';
        Lik(size(options_.varobs,1)+1).isam = 1;
        Lik(size(options_.varobs,1)+1).data = -logpo2;
    end
    Sam.name = bayestopt_.name;
    Sam.dim  = [size(x) 0];
    Sam.data = [x];
    
    Rem.id = 'Original';
    Rem.ind= [1:size(x,1)];
    
    Info.dynare=M_.fname;
    Info.order_var=dr_.order_var;
    Out=Out1;
    if options_gsa.ppost
        %     Info.dynare=M_.fname;
        %     Info.order_var=dr_.order_var;
        %     Out=Out1;
        Info.TypeofSample='post';
        save([OutputDirectoryName,'/',fname_,'_glue_post.mat'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info', 'Exo')
        %save([fname_,'_post_glue_smooth'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info')
        
    else
        if options_gsa.pprior
            Info.TypeofSample='prior';
            save([OutputDirectoryName,'/',fname_,'_glue_prior.mat'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info', 'Exo')
            %       save([OutputDirectoryName,'/',fname_,'_prior_glue'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem')
            %       Out=Out1;
            %       save([OutputDirectoryName,'/',fname_,'_prior_glue_smooth'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem')
        else
            Info.TypeofSample='mc';
            save([OutputDirectoryName,'/',fname_,'_glue_mc.mat'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem','Info', 'Exo')
            %       save([OutputDirectoryName,'/',fname_,'_mc_glue'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem')
            %       Out=Out1;
            %       save([OutputDirectoryName,'/',fname_,'_mc_glue_smooth'], 'Out', 'Sam', 'Lik', 'Obs', 'Rem')
        end
    end
    
end
