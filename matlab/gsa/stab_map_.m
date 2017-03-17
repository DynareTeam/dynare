function x0 = stab_map_(OutputDirectoryName,opt_gsa)
%
% function x0 = stab_map_(OutputDirectoryName)
%
% Mapping of stability regions in the prior ranges applying
% Monte Carlo filtering techniques.
%
% INPUTS (from opt_gsa structure)
% Nsam = MC sample size
% fload = 0 to run new MC; 1 to load prevoiusly generated analysis
% alpha2 =  significance level for bivariate sensitivity analysis
% [abs(corrcoef) > alpha2]
% prepSA = 1: save transition matrices for mapping reduced form
%        = 0: no transition matrix saved (default)
% pprior = 1: sample from prior ranges (default): sample saved in
%            _prior.mat   file
%        = 0: sample from posterior ranges: sample saved in
%            _mc.mat file
% OUTPUT:
% x0: one parameter vector for which the model is stable.
%
% GRAPHS
% 1) Pdf's of marginal distributions under the stability (dotted
%     lines) and unstability (solid lines) regions
% 2) Cumulative distributions of:
%   - stable subset (dotted lines)
%   - unacceptable subset (solid lines)
% 3) Bivariate plots of significant correlation patterns
%  ( abs(corrcoef) > alpha2) under the stable and unacceptable subsets
%
% USES qmc_sequence, stab_map_1, stab_map_2
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

%global bayestopt_ estim_params_ dr_ options_ ys_ fname_
global bayestopt_ estim_params_ options_ oo_ M_

% opt_gsa=options_.opt_gsa;

Nsam   = opt_gsa.Nsam;
fload  = opt_gsa.load_stab;
ksstat = opt_gsa.ksstat;
alpha2 = opt_gsa.alpha2_stab;
pvalue_ks = opt_gsa.pvalue_ks;
pvalue_corr = opt_gsa.pvalue_corr;
prepSA = (opt_gsa.redform | opt_gsa.identification);
pprior = opt_gsa.pprior;
neighborhood_width = opt_gsa.neighborhood_width;
ilptau = opt_gsa.ilptau;
nliv   = opt_gsa.morris_nliv;
ntra   = opt_gsa.morris_ntra;

dr_ = oo_.dr;
%if isfield(dr_,'ghx'),
ys_ = oo_.dr.ys;
nspred = M_.nspred; %size(dr_.ghx,2);
nboth = M_.nboth;
nfwrd = M_.nfwrd;
%end
fname_ = M_.fname;

np = estim_params_.np;
nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;
lpmat0=zeros(Nsam,0);
xparam1=[];

pshape = bayestopt_.pshape(nshock+1:end);
p1 = bayestopt_.p1(nshock+1:end);
p2 = bayestopt_.p2(nshock+1:end);
p3 = bayestopt_.p3(nshock+1:end);
p4 = bayestopt_.p4(nshock+1:end);

[junk1,junk2,junk3,lb,ub,junk4] = set_prior(estim_params_,M_,options_); %Prepare bounds
if ~isempty(bayestopt_) && any(bayestopt_.pshape > 0)
    % Set prior bounds
    bounds = prior_bounds(bayestopt_, options_.prior_trunc);
    bounds.lb = max(bounds.lb,lb);
    bounds.ub = min(bounds.ub,ub);
else  % estimated parameters but no declared priors
    % No priors are declared so Dynare will estimate the model by
    % maximum likelihood with inequality constraints for the parameters.
    bounds.lb = lb;
    bounds.ub = ub;
    if opt_gsa.prior_range==0
        warning('GSA:: When using ML, sampling from the prior is not possible. Setting prior_range=1')
        opt_gsa.prior_range=1;
    end
end

if nargin==0,
    OutputDirectoryName='';
end

options_mcf.pvalue_ks = pvalue_ks;
options_mcf.pvalue_corr = pvalue_corr;
options_mcf.alpha2 = alpha2;

name=cell(np,1);
name_tex=cell(np,1);
for jj=1:np        
    if options_.TeX
        [param_name_temp, param_name_tex_temp]= get_the_name(nshock+jj,options_.TeX,M_,estim_params_,options_);
        name_tex{jj,1} = strrep(param_name_tex_temp,'$','');
        name{jj,1} = param_name_temp;
    else
        param_name_temp = get_the_name(nshock+jj,options_.TeX,M_,estim_params_,options_);
        name{jj,1} = param_name_temp;
    end
end
if options_.TeX
    options_mcf.param_names_tex=char(name_tex);
end
options_mcf.param_names = char(name);

options_mcf.fname_ = fname_;
options_mcf.OutputDirectoryName = OutputDirectoryName;
options_mcf.xparam1 = [];

opt=options_;
options_.periods=0;
options_.nomoments=1;
options_.irf=0;
options_.noprint=1;
if fload==0,
    %   if prepSA
    %     T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),Nsam/2);
    %   end

    if isfield(dr_,'ghx'),
        egg=zeros(length(dr_.eigval),Nsam);
    end
    yys=zeros(length(dr_.ys),Nsam);

    if opt_gsa.morris == 1
        [lpmat, OutFact] = Sampling_Function_2(nliv, np+nshock, ntra, ones(np+nshock, 1), zeros(np+nshock,1), []);
        lpmat = lpmat.*(nliv-1)/nliv+1/nliv/2;
        Nsam=size(lpmat,1);
        lpmat0 = lpmat(:,1:nshock);
        lpmat = lpmat(:,nshock+1:end);
%     elseif opt_gsa.morris==3,
%         lpmat = prep_ide(Nsam,np,5);
%         Nsam=size(lpmat,1);
    else
        if np<52 && ilptau>0,
            [lpmat] = qmc_sequence(np, int64(1), 0, Nsam)';
            if np>30 || ilptau==2, % scrambled lptau
                for j=1:np,
                    lpmat(:,j)=lpmat(randperm(Nsam),j);
                end
            end
        else %ilptau==0
            [lpmat] = NaN(Nsam,np);
            for j=1:np,
                lpmat(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
            end

        end
    end
    %   try
    dummy=prior_draw_gsa(1); %initialize persistent variables
    %   catch
    %     if pprior,
    %       if opt_gsa.prior_range==0;
    %         error('Some unknown prior is specified or ML estimation,: use prior_range=1 option!!');
    %       end
    %     end
    %
    %   end
    if pprior,
        for j=1:nshock,
            if opt_gsa.morris~=1,
                lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
            end
            if opt_gsa.prior_range
                lpmat0(:,j)=lpmat0(:,j).*(bounds.ub(j)-bounds.lb(j))+bounds.lb(j);
            end
        end
        if opt_gsa.prior_range
            %       if opt_gsa.identification,
            %         deltx=min(0.001, 1/Nsam/2);
            %         for j=1:np,
            %           xdelt(:,:,j)=prior_draw_gsa(0,[lpmat0 lpmat]+deltx);
            %         end
            %       end
            for j=1:np,
                lpmat(:,j)=lpmat(:,j).*(bounds.ub(j+nshock)-bounds.lb(j+nshock))+bounds.lb(j+nshock);
            end
        else
            xx=prior_draw_gsa(0,[lpmat0 lpmat]);
            %       if opt_gsa.identification,
            %         deltx=min(0.001, 1/Nsam/2);
            %         ldum=[lpmat0 lpmat];
            %         ldum = prior_draw_gsa(0,ldum+deltx);
            %         for j=1:nshock+np,
            %           xdelt(:,:,j)=xx;
            %           xdelt(:,j,j)=ldum(:,j);
            %         end
            %         clear ldum
            %       end
            lpmat0=xx(:,1:nshock);
            lpmat=xx(:,nshock+1:end);
            clear xx;
        end
    else
        %         for j=1:nshock,
        %             xparam1(j) = oo_.posterior_mode.shocks_std.(bayestopt_.name{j});
        %             sd(j) = oo_.posterior_std.shocks_std.(bayestopt_.name{j});
        %             lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
        %             lb = max(bayestopt_.lb(j), xparam1(j)-2*sd(j));
        %             ub1=xparam1(j)+(xparam1(j) - lb); % define symmetric range around the mode!
        %             ub = min(bayestopt_.ub(j),ub1);
        %             if ub<ub1,
        %                 lb=xparam1(j)-(ub-xparam1(j)); % define symmetric range around the mode!
        %             end
        %             lpmat0(:,j) = lpmat0(:,j).*(ub-lb)+lb;
        %         end
        %         %
        %         for j=1:np,
        %             xparam1(j+nshock) = oo_.posterior_mode.parameters.(bayestopt_.name{j+nshock});
        %             sd(j+nshock) = oo_.posterior_std.parameters.(bayestopt_.name{j+nshock});
        %             lb = max(bayestopt_.lb(j+nshock),xparam1(j+nshock)-2*sd(j+nshock));
        %             ub1=xparam1(j+nshock)+(xparam1(j+nshock) - lb); % define symmetric range around the mode!
        %             ub = min(bayestopt_.ub(j+nshock),ub1);
        %             if ub<ub1,
        %                 lb=xparam1(j+nshock)-(ub-xparam1(j+nshock)); % define symmetric range around the mode!
        %             end
        %             %ub = min(bayestopt_.ub(j+nshock),xparam1(j+nshock)+2*sd(j+nshock));
        %             if np>30 & np<52
        %                 lpmat(:,j) = lpmat(randperm(Nsam),j).*(ub-lb)+lb;
        %             else
        %                 lpmat(:,j) = lpmat(:,j).*(ub-lb)+lb;
        %             end
        %         end
        %load([fname_,'_mode'])
        if neighborhood_width>0 && isempty(options_.mode_file),
            xparam1 = get_all_parameters(estim_params_,M_);
        else
        eval(['load ' options_.mode_file '.mat;']);
        end
        if neighborhood_width>0,
            for j=1:nshock,
                if opt_gsa.morris ~= 1
                   lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
                end
                ub=min([bounds.ub(j) xparam1(j)*(1+neighborhood_width)]);
                lb=max([bounds.lb(j) xparam1(j)*(1-neighborhood_width)]);
                lpmat0(:,j)=lpmat0(:,j).*(ub-lb)+lb;
            end
            for j=1:np,
                ub=xparam1(j+nshock)*(1+sign(xparam1(j+nshock))*neighborhood_width);
                lb=xparam1(j+nshock)*(1-sign(xparam1(j+nshock))*neighborhood_width);
                if bounds.ub(j+nshock)>=xparam1(j) && bounds.lb(j)<=xparam1(j+nshock),
                    ub=min([bounds.ub(j+nshock) ub]);
                    lb=max([bounds.lb(j+nshock) lb]);
                else
                    fprintf('\nstab_map_:: the calibrated value of param %s for neighborhood_width sampling is outside prior bounds.\nWe allow violation of bounds for this parameter, but if this was not done on purpose, please change calibration before running neighborhood_width sampling\n', bayestopt_.name{j+nshock})
                end
                lpmat(:,j)=lpmat(:,j).*(ub-lb)+lb;
            end
        else
            d = chol(inv(hh));
            lp=randn(Nsam*2,nshock+np)*d+kron(ones(Nsam*2,1),xparam1');
            lnprior=zeros(1,Nsam*2);
            for j=1:Nsam*2,
                lnprior(j) = any(lp(j,:)'<=bounds.lb | lp(j,:)'>=bounds.ub);
            end
            ireal=[1:2*Nsam];
            ireal=ireal(find(lnprior==0));
            lp=lp(ireal,:);
            Nsam=min(Nsam, length(ireal));
            lpmat0=lp(1:Nsam,1:nshock);
            lpmat=lp(1:Nsam,nshock+1:end);
            clear lp lnprior ireal;
        end
    end
    %
    h = dyn_waitbar(0,'Please wait...');
    istable=[1:Nsam];
    jstab=0;
    iunstable=[1:Nsam];
    iindeterm=zeros(1,Nsam);
    iwrong=zeros(1,Nsam);
    inorestriction=zeros(1,Nsam);
    irestriction=zeros(1,Nsam);
    infox=zeros(Nsam,1);
    for j=1:Nsam,
        M_ = set_all_parameters([lpmat0(j,:) lpmat(j,:)]',estim_params_,M_);
        %try stoch_simul([]);
        try
            if ~ isempty(options_.endogenous_prior_restrictions.moment)
                [Tt,Rr,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);
            else
                [Tt,Rr,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_,'restrict');
            end
            infox(j,1)=info(1);
            if infox(j,1)==0 && ~exist('T','var'),
                dr_=oo_.dr;
                if prepSA,
                    try
                        T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),Nsam);
                    catch
                        ME = lasterror();
                        if strcmp('MATLAB:nomem',ME.identifier),
                            prepSA=0;
                            disp('The model is too large for storing state space matrices ...')
                            disp('for mapping reduced form or for identification')
                        end
                        T=[];
                    end
                else
                    T=[];
                end
                egg=zeros(length(dr_.eigval),Nsam);
            end
            if infox(j,1),
%                 disp('no solution'),
                if isfield(oo_.dr,'ghx'),
                    oo_.dr=rmfield(oo_.dr,'ghx');
                end
                if (infox(j,1)<3 || infox(j,1)>5) && isfield(oo_.dr,'eigval'),
                    oo_.dr=rmfield(oo_.dr,'eigval');
                end
            end
        catch ME
            if isfield(oo_.dr,'eigval'),
                oo_.dr=rmfield(oo_.dr,'eigval');
            end
            if isfield(oo_.dr,'ghx'),
                oo_.dr=rmfield(oo_.dr,'ghx');
            end
            disp('No solution could be found'),
        end
        dr_ = oo_.dr;
        if isfield(dr_,'ghx'),
            egg(:,j) = sort(dr_.eigval);
            if prepSA
                jstab=jstab+1;
                T(:,:,jstab) = [dr_.ghx dr_.ghu];
                %         [A,B] = ghx2transition(squeeze(T(:,:,jstab)), ...
                %           bayestopt_.restrict_var_list, ...
                %           bayestopt_.restrict_columns, ...
                %           bayestopt_.restrict_aux);
            end
            if ~exist('nspred','var'),
                nspred = dr_.nspred; %size(dr_.ghx,2);
                nboth = dr_.nboth;
                nfwrd = dr_.nfwrd;
            end
            info=endogenous_prior_restrictions(Tt,Rr,M_,options_,oo_);
            infox(j,1)=info(1);
            if info(1),
                inorestriction(j)=j;
            else
                iunstable(j)=0;
                irestriction(j)=j;
            end
        else
            istable(j)=0;
            if isfield(dr_,'eigval')
                egg(:,j) = sort(dr_.eigval);
                if exist('nspred','var')
                    if any(isnan(egg(1:nspred,j)))
                        iwrong(j)=j;
                    else
                        if (nboth || nfwrd) && abs(egg(nspred+1,j))<=options_.qz_criterium,
                            iindeterm(j)=j;
                        end
                    end
                end
            else
                if exist('egg','var'),
                    egg(:,j)=ones(size(egg,1),1).*NaN;
                end
                iwrong(j)=j;
            end
        end
        ys_=real(dr_.ys);
        yys(:,j) = ys_;
        ys_=yys(:,1);
        dyn_waitbar(j/Nsam,h,['MC iteration ',int2str(j),'/',int2str(Nsam)])
    end
    dyn_waitbar_close(h);
    if prepSA && jstab,
        T=T(:,:,1:jstab);
    else
        T=[];
    end
    istable=istable(find(istable));  % stable params ignoring restrictions
    irestriction=irestriction(find(irestriction));  % stable params & restrictions OK
    inorestriction=inorestriction(find(inorestriction));  % stable params violating restrictions
    iunstable=iunstable(find(iunstable));   % violation of BK & restrictions & solution could not be found (whatever goes wrong)
    iindeterm=iindeterm(find(iindeterm));  % indeterminacy
    iwrong=iwrong(find(iwrong));  % dynare could not find solution
    ixun=iunstable(find(~ismember(iunstable,[iindeterm,iwrong,inorestriction]))); % explosive roots

    %     % map stable samples
    %     istable=[1:Nsam];
    %     for j=1:Nsam,
    %         if any(isnan(egg(1:nspred,j)))
    %             istable(j)=0;
    %         else
    %             if abs(egg(nspred,j))>=options_.qz_criterium; %(1-(options_.qz_criterium-1)); %1-1.e-5;
    %                 istable(j)=0;
    %                 %elseif (dr_.nboth | dr_.nfwrd) & abs(egg(nspred+1,j))<=options_.qz_criterium; %1+1.e-5;
    %             elseif (nboth | nfwrd) & abs(egg(nspred+1,j))<=options_.qz_criterium; %1+1.e-5;
    %                 istable(j)=0;
    %             end
    %         end
    %     end
    %     istable=istable(find(istable));  % stable params
    %
    %     % map unstable samples
    %     iunstable=[1:Nsam];
    %     for j=1:Nsam,
    %         %if abs(egg(dr_.npred+1,j))>1+1.e-5 & abs(egg(dr_.npred,j))<1-1.e-5;
    %         %if (dr_.nboth | dr_.nfwrd),
    %         if ~any(isnan(egg(1:5,j)))
    %             if (nboth | nfwrd),
    %                 if abs(egg(nspred+1,j))>options_.qz_criterium & abs(egg(nspred,j))<options_.qz_criterium; %(1-(options_.qz_criterium-1));
    %                     iunstable(j)=0;
    %                 end
    %             else
    %                 if abs(egg(nspred,j))<options_.qz_criterium; %(1-(options_.qz_criterium-1));
    %                     iunstable(j)=0;
    %                 end
    %             end
    %         end
    %     end
    %     iunstable=iunstable(find(iunstable));   % unstable params
    bkpprior.pshape=bayestopt_.pshape;
    bkpprior.p1=bayestopt_.p1;
    bkpprior.p2=bayestopt_.p2;
    bkpprior.p3=bayestopt_.p3;
    bkpprior.p4=bayestopt_.p4;
    if pprior,
        if ~prepSA
            save([OutputDirectoryName filesep fname_ '_prior.mat'], ...
                'bkpprior','lpmat','lpmat0','irestriction','iunstable','istable','iindeterm','iwrong','ixun', ...
                'egg','yys','nspred','nboth','nfwrd','infox')
        else
            save([OutputDirectoryName filesep fname_ '_prior.mat'], ...
                'bkpprior','lpmat','lpmat0','irestriction','iunstable','istable','iindeterm','iwrong','ixun', ...
                'egg','yys','T','nspred','nboth','nfwrd','infox')
        end

    else
        if ~prepSA
            save([OutputDirectoryName filesep fname_ '_mc.mat'], ...
                'lpmat','lpmat0','irestriction','iunstable','istable','iindeterm','iwrong','ixun', ...
                'egg','yys','nspred','nboth','nfwrd','infox')
        else
            save([OutputDirectoryName filesep fname_ '_mc.mat'], ...
                'lpmat','lpmat0','irestriction','iunstable','istable','iindeterm','iwrong','ixun', ...
                'egg','yys','T','nspred','nboth','nfwrd','infox')
        end
    end
else
    if pprior,
        filetoload=[OutputDirectoryName filesep fname_ '_prior.mat'];
    else
        filetoload=[OutputDirectoryName filesep fname_ '_mc.mat'];
    end
    load(filetoload,'lpmat','lpmat0','irestriction','iunstable','istable','iindeterm','iwrong','ixun','egg','yys','nspred','nboth','nfwrd','infox')
    Nsam = size(lpmat,1);
    if pprior==0 && ~isempty(options_.mode_file),
        eval(['load ' options_.mode_file '.mat;']);
    end


    if prepSA && isempty(strmatch('T',who('-file', filetoload),'exact')),
        h = dyn_waitbar(0,'Please wait...');
        options_.periods=0;
        options_.nomoments=1;
        options_.irf=0;
        options_.noprint=1;
        stoch_simul([]);
        %T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),length(istable));
        ntrans=length(istable);
        yys=NaN(length(ys_),ntrans);
        for j=1:ntrans,
            M_.params(estim_params_.param_vals(:,1)) = lpmat(istable(j),:)';
            %stoch_simul([]);
            [Tt,Rr,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_,'restrict');
            % This syntax is not compatible with the current version of dynare_resolve [stepan].
            %[Tt,Rr,SteadyState,info] = dynare_resolve(bayestopt_.restrict_var_list,...
            %    bayestopt_.restrict_columns,...
            %    bayestopt_.restrict_aux);
            if ~exist('T','var')
                T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),ntrans);
            end
            dr_ = oo_.dr;
            T(:,:,j) = [dr_.ghx dr_.ghu];
            if ~exist('nspred','var')
                nspred = dr_.nspred; %size(dr_.ghx,2);
                nboth = dr_.nboth;
                nfwrd = dr_.nfwrd;
            end
            ys_=real(dr_.ys);
            yys(:,j) = ys_;
            ys_=yys(:,1);
            dyn_waitbar(j/ntrans,h,['MC iteration ',int2str(j),'/',int2str(ntrans)])
        end
        dyn_waitbar_close(h);
        save(filetoload,'T','-append')
    elseif prepSA
        load(filetoload,'T')
    end
end

if pprior
    aunstname='prior_unstable'; aunsttitle='Prior StabMap: explosiveness of solution';
    aindname='prior_indeterm'; aindtitle='Prior StabMap: Indeterminacy';
    awrongname='prior_wrong'; awrongtitle='Prior StabMap: inability to find solution';
    acalibname='prior_calib'; acalibtitle='Prior StabMap: IRF/moment restrictions';
    asname='prior_stable'; atitle='Prior StabMap: Parameter driving non-existence of unique stable solution (Unacceptable)';
else
    aunstname='mc_unstable'; aunsttitle='MC (around posterior mode) StabMap: explosiveness of solution';
    aindname='mc_indeterm';  aindtitle='MC (around posterior mode) StabMap: Indeterminacy';
    awrongname='mc_wrong'; awrongtitle='MC (around posterior mode) StabMap: inability to find solution';
    acalibname='mc_calib'; acalibtitle='MC (around posterior mode) StabMap: IRF/moment restrictions';
    asname='mc_stable'; atitle='MC (around posterior mode) StabMap: Parameter driving non-existence of unique stable solution (Unacceptable)';
end
delete([OutputDirectoryName,filesep,fname_,'_',asname,'.*']);
delete([OutputDirectoryName,filesep,fname_,'_',acalibname,'.*']);
delete([OutputDirectoryName,filesep,fname_,'_',aindname,'.*']);
delete([OutputDirectoryName,filesep,fname_,'_',aunstname,'.*']);
delete([OutputDirectoryName,filesep,fname_,'_',awrongname,'.*']);

if length(iunstable)>0 || length(iwrong)>0,
    fprintf(['%4.1f%% of the prior support gives unique saddle-path solution.\n'],length(istable)/Nsam*100)
    fprintf(['%4.1f%% of the prior support gives explosive dynamics.\n'],(length(ixun) )/Nsam*100)
    if ~isempty(iindeterm),
        fprintf(['%4.1f%% of the prior support gives indeterminacy.\n'],length(iindeterm)/Nsam*100)
    end
    inorestriction = istable(find(~ismember(istable,irestriction))); % violation of prior restrictions
    if ~isempty(iwrong) || ~isempty(inorestriction),
        skipline()
        if any(infox==49),
            fprintf(['%4.1f%% of the prior support violates prior restrictions.\n'],(length(inorestriction) )/Nsam*100)
        end
        if ~isempty(iwrong),
            skipline()
            disp(['For ',num2str(length(iwrong)/Nsam*100,'%4.1f'),'% of the prior support dynare could not find a solution.'])
            skipline()
        end
        if any(infox==1),
            disp(['    For ',num2str(length(find(infox==1))/Nsam*100,'%4.1f'),'% The model doesn''t determine the current variables uniquely.'])
        end
        if any(infox==2),
            disp(['    For ',num2str(length(find(infox==2))/Nsam*100,'%4.1f'),'% MJDGGES returned an error code.'])
        end
        if any(infox==6),
            disp(['    For ',num2str(length(find(infox==6))/Nsam*100,'%4.1f'),'% The jacobian evaluated at the deterministic steady state is complex.'])
        end
        if any(infox==19),
            disp(['    For ',num2str(length(find(infox==19))/Nsam*100,'%4.1f'),'% The steadystate routine thrown an exception (inconsistent deep parameters).'])
        end
        if any(infox==20),
            disp(['    For ',num2str(length(find(infox==20))/Nsam*100,'%4.1f'),'% Cannot find the steady state.'])
        end
        if any(infox==21),
            disp(['    For ',num2str(length(find(infox==21))/Nsam*100,'%4.1f'),'% The steady state is complex.'])
        end
        if any(infox==22),
            disp(['    For ',num2str(length(find(infox==22))/Nsam*100,'%4.1f'),'% The steady has NaNs.'])
        end
        if any(infox==23),
            disp(['    For ',num2str(length(find(infox==23))/Nsam*100,'%4.1f'),'% M_.params has been updated in the steadystate routine and has complex valued scalars.'])
        end
        if any(infox==24),
            disp(['    For ',num2str(length(find(infox==24))/Nsam*100,'%4.1f'),'% M_.params has been updated in the steadystate routine and has some NaNs.'])
        end
        if any(infox==30),
            disp(['    For ',num2str(length(find(infox==30))/Nsam*100,'%4.1f'),'% Ergodic variance can''t be computed.'])
        end

    end
    skipline()
    if length(iunstable)<Nsam || length(istable)>1
        itot = [1:Nsam];
        isolve = itot(find(~ismember(itot,iwrong))); % dynare could find a solution
        % Blanchard Kahn
        if neighborhood_width,
            options_mcf.xparam1 = xparam1(nshock+1:end);
        end
        itmp = itot(find(~ismember(itot,istable)));
        options_mcf.amcf_name = asname;
        options_mcf.amcf_title = atitle;
        options_mcf.beha_title = 'unique Stable Saddle-Path';
        options_mcf.nobeha_title = 'NO unique Stable Saddle-Path';
        options_mcf.title = 'unique solution';
        mcf_analysis(lpmat, istable, itmp, options_mcf, options_);

        if ~isempty(iindeterm),
            itmp = isolve(find(~ismember(isolve,iindeterm)));
            options_mcf.amcf_name = aindname;
            options_mcf.amcf_title = aindtitle;
            options_mcf.beha_title = 'NO indeterminacy';
            options_mcf.nobeha_title = 'indeterminacy';
            options_mcf.title = 'indeterminacy';
            mcf_analysis(lpmat, itmp, iindeterm, options_mcf, options_);
        end
        
        if ~isempty(ixun),
            itmp = isolve(find(~ismember(isolve,ixun)));
            options_mcf.amcf_name = aunstname;
            options_mcf.amcf_title = aunsttitle;
            options_mcf.beha_title = 'NO explosive solution';
            options_mcf.nobeha_title = 'explosive solution';
            options_mcf.title = 'instability';
            mcf_analysis(lpmat, itmp, ixun, options_mcf, options_);
        end
        
        inorestriction = istable(find(~ismember(istable,irestriction))); % violation of prior restrictions
        iwrong = iwrong(find(~ismember(iwrong,inorestriction))); % what went wrong beyond prior restrictions
        if ~isempty(iwrong),
            itmp = itot(find(~ismember(itot,iwrong)));
            options_mcf.amcf_name = awrongname;
            options_mcf.amcf_title = awrongtitle;
            options_mcf.beha_title = 'NO inability to find a solution';
            options_mcf.nobeha_title = 'inability to find a solution';
            options_mcf.title = 'inability to find a solution';
            mcf_analysis(lpmat, itmp, iwrong, options_mcf, options_);
        end
        
        if ~isempty(irestriction),
            if neighborhood_width,
                options_mcf.xparam1 = xparam1;
            end
            np=size(bayestopt_.name,1);
            name=cell(np,1);
            name_tex=cell(np,1);
            for jj=1:np
                if options_.TeX
                    [param_name_temp, param_name_tex_temp]= get_the_name(jj,options_.TeX,M_,estim_params_,options_);
                    name_tex{jj,1} = strrep(param_name_tex_temp,'$','');
                    name{jj,1} = param_name_temp;
                else
                    param_name_temp = get_the_name(jj,options_.TeX,M_,estim_params_,options_);
                    name{jj,1} = param_name_temp;
                end
            end
            if options_.TeX
                options_mcf.param_names_tex = char(name_tex);
            end
            options_mcf.param_names = char(name);
            options_mcf.amcf_name = acalibname;
            options_mcf.amcf_title = acalibtitle;
            options_mcf.beha_title = 'prior IRF/moment calibration';
            options_mcf.nobeha_title = 'NO prior IRF/moment calibration';
            options_mcf.title = 'prior restrictions';
            mcf_analysis([lpmat0 lpmat], irestriction, inorestriction, options_mcf, options_);
            iok = irestriction(1);
            x0 = [lpmat0(iok,:)'; lpmat(iok,:)'];
        else
            iok = istable(1);
            x0=0.5.*(bounds.ub(1:nshock)-bounds.lb(1:nshock))+bounds.lb(1:nshock);
            x0 = [x0; lpmat(iok,:)'];
        end
        
        M_ = set_all_parameters(x0,estim_params_,M_);
        [oo_.dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        %     stoch_simul([]);
    else
        disp('All parameter values in the specified ranges are not acceptable!')
        x0=[];
    end
else
    disp('All parameter values in the specified ranges give unique saddle-path solution,')
    disp('and match prior IRF/moment restriction(s) if any!')
    x0=0.5.*(bounds.ub(1:nshock)-bounds.lb(1:nshock))+bounds.lb(1:nshock);
    x0 = [x0; lpmat(istable(1),:)'];

end

xparam1=x0;
save prior_ok.mat xparam1;

options_.periods=opt.periods;
if isfield(opt,'nomoments'),
    options_.nomoments=opt.nomoments;
end
options_.irf=opt.irf;
options_.noprint=opt.noprint;




