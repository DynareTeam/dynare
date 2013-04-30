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
lpmat0=[];
xparam1=[];

pshape = bayestopt_.pshape(nshock+1:end);
p1 = bayestopt_.p1(nshock+1:end);
p2 = bayestopt_.p2(nshock+1:end);
p3 = bayestopt_.p3(nshock+1:end);
p4 = bayestopt_.p4(nshock+1:end);

if nargin==0,
    OutputDirectoryName='';
end

opt=options_;
options_.periods=0;
options_.nomoments=1;
options_.irf=0;
options_.noprint=1;
options_.simul=0;
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
            %[lpmat] = rand(Nsam,np);
            for j=1:np,
                lpmat(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
            end

        end
    end
    %   try
    dummy=prior_draw_gsa(1);
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
                lpmat0(:,j)=lpmat0(:,j).*(bayestopt_.ub(j)-bayestopt_.lb(j))+bayestopt_.lb(j);
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
                lpmat(:,j)=lpmat(:,j).*(bayestopt_.ub(j+nshock)-bayestopt_.lb(j+nshock))+bayestopt_.lb(j+nshock);
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
        eval(['load ' options_.mode_file '.mat;']);
        if neighborhood_width>0,
            for j=1:nshock,
                lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
                ub=min([bayestopt_.ub(j) xparam1(j)*(1+neighborhood_width)]);
                lb=max([bayestopt_.lb(j) xparam1(j)*(1-neighborhood_width)]);
                lpmat0(:,j)=lpmat0(:,j).*(ub-lb)+lb;
            end
            for j=1:np,
                ub=min([bayestopt_.ub(j+nshock) xparam1(j+nshock)*(1+neighborhood_width)]);
                lb=max([bayestopt_.lb(j+nshock) xparam1(j+nshock)*(1-neighborhood_width)]);
                lpmat(:,j)=lpmat(:,j).*(ub-lb)+lb;
            end
        else
            d = chol(inv(hh));
            lp=randn(Nsam*2,nshock+np)*d+kron(ones(Nsam*2,1),xparam1');
            for j=1:Nsam*2,
                lnprior(j) = any(lp(j,:)'<=bayestopt_.lb | lp(j,:)'>=bayestopt_.ub);
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
    for j=1:Nsam,
        M_.params(estim_params_.param_vals(:,1)) = lpmat(j,:)';
        %try stoch_simul([]);
        try
            [Tt,Rr,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_,'restrict');
            infox(j,1)=info(1);
            if infox(j,1)==0 && ~exist('T'),
                dr_=oo_.dr;
                if prepSA,
                    try
                        T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),Nsam);
                    catch ME
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
            iunstable(j)=0;
            if prepSA
                jstab=jstab+1;
                T(:,:,jstab) = [dr_.ghx dr_.ghu];
                %         [A,B] = ghx2transition(squeeze(T(:,:,jstab)), ...
                %           bayestopt_.restrict_var_list, ...
                %           bayestopt_.restrict_columns, ...
                %           bayestopt_.restrict_aux);
            end
            if ~exist('nspred'),
                nspred = dr_.nspred; %size(dr_.ghx,2);
                nboth = dr_.nboth;
                nfwrd = dr_.nfwrd;
            end
        else
            istable(j)=0;
            if isfield(dr_,'eigval')
                egg(:,j) = sort(dr_.eigval);
                if exist('nspred')
                    if any(isnan(egg(1:nspred,j)))
                        iwrong(j)=j;
                    else
                        if (nboth || nfwrd) && abs(egg(nspred+1,j))<=options_.qz_criterium,
                            iindeterm(j)=j;
                        end
                    end
                end
            else
                if exist('egg'),
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
    istable=istable(find(istable));  % stable params
    iunstable=iunstable(find(iunstable));   % unstable params
    iindeterm=iindeterm(find(iindeterm));  % indeterminacy
    iwrong=iwrong(find(iwrong));  % dynare could not find solution

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
                'bkpprior','lpmat','lpmat0','iunstable','istable','iindeterm','iwrong', ...
                'egg','yys','nspred','nboth','nfwrd','infox')
        else
            save([OutputDirectoryName filesep fname_ '_prior.mat'], ...
                'bkpprior','lpmat','lpmat0','iunstable','istable','iindeterm','iwrong', ...
                'egg','yys','T','nspred','nboth','nfwrd','infox')
        end

    else
        if ~prepSA
            save([OutputDirectoryName filesep fname_ '_mc.mat'], ...
                'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong', ...
                'egg','yys','nspred','nboth','nfwrd','infox')
        else
            save([OutputDirectoryName filesep fname_ '_mc.mat'], ...
                'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong', ...
                'egg','yys','T','nspred','nboth','nfwrd','infox')
        end
    end
else
    if pprior,
        filetoload=[OutputDirectoryName filesep fname_ '_prior.mat'];
    else
        filetoload=[OutputDirectoryName filesep fname_ '_mc.mat'];
    end
    load(filetoload,'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','nspred','nboth','nfwrd','infox')
    Nsam = size(lpmat,1);
    if pprior==0,
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
        for j=1:ntrans,
            M_.params(estim_params_.param_vals(:,1)) = lpmat(istable(j),:)';
            %stoch_simul([]);
            [Tt,Rr,SteadyState,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_,'restrict');
            % This syntax is not compatible with the current version of dynare_resolve [stepan].
            %[Tt,Rr,SteadyState,info] = dynare_resolve(bayestopt_.restrict_var_list,...
            %    bayestopt_.restrict_columns,...
            %    bayestopt_.restrict_aux);
            if ~exist('T')
                T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),ntrans);
            end
            dr_ = oo_.dr;
            T(:,:,j) = [dr_.ghx dr_.ghu];
            if ~exist('nspred')
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
    % univariate
    aname='prior_stab'; atitle='Prior StabMap: Parameter driving non-existence of unique stable solution (Unacceptable)';
    aindetname=[aname, '_indet']; aindettitle='Prior StabMap: Parameter driving indeterminacy';
    aunstablename=[aname, '_unst'];  aunstabletitle='Prior StabMap: Parameter driving explosiveness of solution';
    awronguniname=[aname, '_wrong']; awrongunititle='Prior StabMap: Parameter driving inability to find solution';
    % bivariate
    auname='prior_unacceptable'; autitle='Prior Unacceptable';
    aunstname='prior_unstable'; aunsttitle='Prior Unstable';
    aindname='prior_indeterm'; aindtitle='Prior Indeterminacy';
    awrongname='prior_wrong'; awrongtitle='Prior No Solution Found';
    asname='prior_stable'; astitle='Prior Stable';
else
    % univariate
    aname='mc_stab'; atitle='Posterior StabMap: Parameter driving non-existence of unique stable solution (Unacceptable)';
    aindetname=[aname, '_indet']; aindettitle='Posterior StabMap: Parameter driving indeterminacy';
    aunstablename=[aname, '_unst'];  aunstabletitle='Posterior StabMap: Parameter driving explosiveness of solution';
    awronguniname=[aname, '_wrong']; awrongunititle='Posterior StabMap: Parameter driving inability to find solution';
    % bivariate
    auname='mc_unacceptable'; autitle='Posterior Unacceptable';
    aunstname='mc_unstable'; aunsttitle='Posterior Unstable';
    aindname='mc_indeterm';  aindtitle='Posterior Indeterminacy';
    awrongname='mc_wrong'; awrongtitle='Posterior No Solution Found';
    asname='mc_stable'; astitle='Posterior Stable';
end
delete([OutputDirectoryName,filesep,fname_,'_',aname,'_*.*']);
%delete([OutputDirectoryName,filesep,fname_,'_',aname,'_SA_*.*']);
delete([OutputDirectoryName,filesep,fname_,'_',asname,'_corr_*.*']);
delete([OutputDirectoryName,filesep,fname_,'_',auname,'_corr_*.*']);
delete([OutputDirectoryName,filesep,fname_,'_',aunstname,'_corr_*.*']);
delete([OutputDirectoryName,filesep,fname_,'_',aindname,'_corr_*.*']);

if length(iunstable)>0 && length(iunstable)<Nsam,
    fprintf(['%4.1f%% of the prior support gives unique saddle-path solution.\n'],length(istable)/Nsam*100)
    fprintf(['%4.1f%% of the prior support gives explosive dynamics.\n'],(length(iunstable)-length(iwrong)-length(iindeterm) )/Nsam*100)
    if ~isempty(iindeterm),
        fprintf(['%4.1f%% of the prior support gives indeterminacy.'],length(iindeterm)/Nsam*100)
    end
    if ~isempty(iwrong),
        disp(' ');
        disp(['For ',num2str(length(iwrong)/Nsam*100,'%1.3f'),'\% of the prior support dynare could not find a solution.'])
        disp(' ');
        if any(infox==1),
            disp(['    For ',num2str(length(find(infox==1))/Nsam*100,'%1.3f'),'\% The model doesn''t determine the current variables uniquely.'])
        end
        if any(infox==2),
            disp(['    For ',num2str(length(find(infox==2))/Nsam*100,'%1.3f'),'\% MJDGGES returned an error code.'])
        end
        if any(infox==6),
            disp(['    For ',num2str(length(find(infox==6))/Nsam*100,'%1.3f'),'\% The jacobian evaluated at the deterministic steady state is complex.'])
        end
        if any(infox==19),
            disp(['    For ',num2str(length(find(infox==19))/Nsam*100,'%1.3f'),'\% The steadystate routine thrown an exception (inconsistent deep parameters).'])
        end
        if any(infox==20),
            disp(['    For ',num2str(length(find(infox==20))/Nsam*100,'%1.3f'),'\% Cannot find the steady state.'])
        end
        if any(infox==21),
            disp(['    For ',num2str(length(find(infox==21))/Nsam*100,'%1.3f'),'\% The steady state is complex.'])
        end
        if any(infox==22),
            disp(['    For ',num2str(length(find(infox==22))/Nsam*100,'%1.3f'),'\% The steady has NaNs.'])
        end
        if any(infox==23),
            disp(['    For ',num2str(length(find(infox==23))/Nsam*100,'%1.3f'),'\% M_.params has been updated in the steadystate routine and has complex valued scalars.'])
        end
        if any(infox==24),
            disp(['    For ',num2str(length(find(infox==24))/Nsam*100,'%1.3f'),'\% M_.params has been updated in the steadystate routine and has some NaNs.'])
        end
        if any(infox==30),
            disp(['    For ',num2str(length(find(infox==30))/Nsam*100,'%1.3f'),'\% Ergodic variance can''t be computed.'])
        end

    end
    disp(' ');
    % Blanchard Kahn
    [proba, dproba] = stab_map_1(lpmat, istable, iunstable, aname,0);
%     indstab=find(dproba>ksstat);
    indstab=find(proba<pvalue_ks);
    disp('Smirnov statistics in driving acceptable behaviour')
    for j=1:length(indstab),
        disp([M_.param_names(estim_params_.param_vals(indstab(j),1),:),'   d-stat = ', num2str(dproba(indstab(j)),'%1.3f'),'   p-value = ', num2str(proba(indstab(j)),'%1.3f')])
    end
    disp(' ');
    if ~isempty(indstab)
        stab_map_1(lpmat, istable, iunstable, aname, 1, indstab, OutputDirectoryName,[],atitle);
    end
    ixun=iunstable(find(~ismember(iunstable,[iindeterm,iwrong])));
    if ~isempty(iindeterm),
        [proba, dproba] = stab_map_1(lpmat, [1:Nsam], iindeterm, aindetname ,0);
%         indindet=find(dproba>ksstat);
        indindet=find(proba<pvalue_ks);
        disp('Smirnov statistics in driving indeterminacy')
        for j=1:length(indindet),
            disp([M_.param_names(estim_params_.param_vals(indindet(j),1),:),'   d-stat = ', num2str(dproba(indindet(j)),'%1.3f'),'   p-value = ', num2str(proba(indindet(j)),'%1.3f')])
        end
        disp(' ');
        if ~isempty(indindet)
            stab_map_1(lpmat, [1:Nsam], iindeterm, aindetname, 1, indindet, OutputDirectoryName,[],aindettitle);
        end
    end

    if ~isempty(ixun),
        [proba, dproba] = stab_map_1(lpmat, [1:Nsam], ixun, aunstablename,0);
%         indunst=find(dproba>ksstat);
        indunst=find(proba<pvalue_ks);
        disp('Smirnov statistics in driving instability')
        for j=1:length(indunst),
            disp([M_.param_names(estim_params_.param_vals(indunst(j),1),:),'   d-stat = ', num2str(dproba(indunst(j)),'%1.3f'),'   p-value = ', num2str(proba(indunst(j)),'%1.3f')])
        end
        disp(' ');
        if ~isempty(indunst)
            stab_map_1(lpmat, [1:Nsam], ixun, aunstablename, 1, indunst, OutputDirectoryName,[],aunstabletitle);
        end
    end

    if ~isempty(iwrong),
        [proba, dproba] = stab_map_1(lpmat, [1:Nsam], iwrong, awronguniname,0);
%         indwrong=find(dproba>ksstat);
        indwrong=find(proba<pvalue_ks);
        disp('Smirnov statistics in driving no solution')
        for j=1:length(indwrong),
            disp([M_.param_names(estim_params_.param_vals(indwrong(j),1),:),'   d-stat = ', num2str(dproba(indwrong(j)),'%1.3f'),'   p-value = ', num2str(proba(indwrong(j)),'%1.3f')])
        end
        disp(' ');
        if ~isempty(indwrong)
            stab_map_1(lpmat, [1:Nsam], iwrong, awronguniname, 1, indwrong, OutputDirectoryName,[],awrongunititle);
        end
    end

    disp(' ')
    disp('Starting bivariate analysis:')

    c0=corrcoef(lpmat(istable,:));
    c00=tril(c0,-1);

    stab_map_2(lpmat(istable,:),alpha2, pvalue_corr, asname, OutputDirectoryName,xparam1,astitle);
    if length(iunstable)>10,
        stab_map_2(lpmat(iunstable,:),alpha2, pvalue_corr, auname, OutputDirectoryName,xparam1,autitle);
    end
    if length(iindeterm)>10,
        stab_map_2(lpmat(iindeterm,:),alpha2, pvalue_corr, aindname, OutputDirectoryName,xparam1,aindtitle);
    end
    if length(ixun)>10,
        stab_map_2(lpmat(ixun,:),alpha2, pvalue_corr, aunstname, OutputDirectoryName,xparam1,aunsttitle);
    end
    if length(iwrong)>10,
        stab_map_2(lpmat(iwrong,:),alpha2, pvalue_corr, awrongname, OutputDirectoryName,xparam1,awrongtitle);
    end

    x0=0.5.*(bayestopt_.ub(1:nshock)-bayestopt_.lb(1:nshock))+bayestopt_.lb(1:nshock);
    x0 = [x0; lpmat(istable(1),:)'];
    if istable(end)~=Nsam
        M_.params(estim_params_.param_vals(:,1)) = lpmat(istable(1),:)';
        [oo_.dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        %     stoch_simul([]);
    end
else
    if length(iunstable)==0,
        disp('All parameter values in the specified ranges give unique saddle-path solution!')
        x0=0.5.*(bayestopt_.ub(1:nshock)-bayestopt_.lb(1:nshock))+bayestopt_.lb(1:nshock);
        x0 = [x0; lpmat(istable(1),:)'];
    else
        disp('All parameter values in the specified ranges are not acceptable!')
        x0=[];
    end

end

xparam1=x0;
save prior_ok xparam1;

options_.periods=opt.periods;
if isfield(opt,'nomoments'),
    options_.nomoments=opt.nomoments;
end
options_.irf=opt.irf;
options_.noprint=opt.noprint;
if isfield(opt,'simul'),
    options_.simul=opt.simul;
end



