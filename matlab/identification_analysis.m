function [ide_hess, ide_moments, ide_model, ide_lre, derivatives_info, info] = identification_analysis(params,indx,indexo,options_ident,data_info, prior_exist, name_tex, init)
% function [ide_hess, ide_moments, ide_model, ide_lre, derivatives_info, info] = identification_analysis(params,indx,indexo,options_ident,data_info, prior_exist, name_tex, init)
% given the parameter vector params, wraps all identification analyses
%
% INPUTS
%    o params             [array] parameter values for identification checks
%    o indx               [array] index of estimated parameters
%    o indexo             [array] index of estimated shocks
%    o options_ident      [structure] identification options
%    o data_info          [structure] data info for Kalman Filter
%    o prior_exist        [integer] 
%                           =1 when prior exists and indentification is checked only for estimated params and shocks
%                           =0 when prior is not defined and indentification is checked for all params and shocks
%    o nem_tex            [char] list of tex names
%    o init               [integer] flag  for initialization of persistent vars
%    
% OUTPUTS
%    o ide_hess           [structure] identification results using Asymptotic Hessian
%    o ide_moments        [structure] identification results using theoretical moments
%    o ide_model          [structure] identification results using reduced form solution
%    o ide_lre            [structure] identification results using LRE model
%    o derivatives_info   [structure] info about analytic derivs
%    o info               output from dynare resolve
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

global oo_ M_ options_ bayestopt_ estim_params_
persistent indH indJJ indLRE

nparam=length(params);
np=length(indx);
offset=nparam-np;
M_ = set_all_parameters(params,estim_params_,M_);

nlags = options_ident.ar;
useautocorr = options_ident.useautocorr;
advanced = options_ident.advanced;
replic = options_ident.replic;
periods = options_ident.periods;
max_dim_cova_group = options_ident.max_dim_cova_group;
normalize_jacobians = options_ident.normalize_jacobians;
kron_flag = options_ident.analytic_derivation_mode;
    
[I,J]=find(M_.lead_lag_incidence');

ide_hess = struct();
ide_moments = struct();
ide_model = struct();
ide_lre = struct();
derivatives_info = struct();

[A,B,ys,info,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);
if info(1)==0,
    oo0=oo_;
    tau=[oo_.dr.ys(oo_.dr.order_var); vec(A); dyn_vech(B*M_.Sigma_e*B')];
    yy0=oo_.dr.ys(I);
    [residual, g1 ] = feval([M_.fname,'_dynamic'],yy0, ...
        oo_.exo_steady_state', M_.params, ...
        oo_.dr.ys, 1);
    vg1 = [oo_.dr.ys(oo_.dr.order_var); vec(g1)];

    [JJ, H, gam, gp, dA, dOm, dYss] = getJJ(A, B, M_,oo0,options_,kron_flag,indx,indexo,bayestopt_.mf2,nlags,useautocorr);
    derivatives_info.DT=dA;
    derivatives_info.DOm=dOm;
    derivatives_info.DYss=dYss;
    if init,
        indJJ = (find(max(abs(JJ'),[],1)>1.e-8));
        while length(indJJ)<nparam && nlags<10,
            disp('The number of moments with non-zero derivative is smaller than the number of parameters')
            disp(['Try increasing ar = ', int2str(nlags+1)])           
            nlags=nlags+1;
            [JJ, H, gam, gp, dA, dOm, dYss] = getJJ(A, B, M_,oo0,options_,kron_flag,indx,indexo,bayestopt_.mf2,nlags,useautocorr);
            derivatives_info.DT=dA;
            derivatives_info.DOm=dOm;
            derivatives_info.DYss=dYss;
            evalin('caller',['options_ident.ar=',int2str(nlags),';']);
        end
        if length(indJJ)<nparam,
            disp('The number of moments with non-zero derivative is smaller than the number of parameters')
            disp('up to 10 lags: check your model')           
            disp('Either further increase ar or reduce the list of estimated parameters')           
            error('IDETooManyParams',''),
        end
        indH = (find(max(abs(H'),[],1)>1.e-8));
        indLRE = (find(max(abs(gp'),[],1)>1.e-8));
    end
    TAU(:,1)=tau(indH);
    LRE(:,1)=vg1(indLRE);
    GAM(:,1)=gam(indJJ);
    siJ = (JJ(indJJ,:));
    siH = (H(indH,:));   
    siLRE = (gp(indLRE,:));
    ide_strength_J=NaN(1,nparam);
    ide_strength_J_prior=NaN(1,nparam);
    if init, %~isempty(indok),
        normaliz = abs(params);
        if prior_exist,
            if ~isempty(estim_params_.var_exo),
                normaliz1 = estim_params_.var_exo(:,7); % normalize with prior standard deviation
            else
                normaliz1=[];
            end
            if ~isempty(estim_params_.param_vals),
                normaliz1 = [normaliz1; estim_params_.param_vals(:,7)]'; % normalize with prior standard deviation
            end
            %                         normaliz = max([normaliz; normaliz1]);
            normaliz1(isinf(normaliz1)) = 1;
            
        else
            normaliz1 = NaN(1,nparam);
        end
        try,
            options_.irf = 0;
            options_.noprint = 1;
            options_.order = 1;
            options_.SpectralDensity.trigger = 0;
            options_.periods = data_info.info.ntobs+100;
            if options_.kalman_algo > 2,
                options_.kalman_algo = 1;
            end
            analytic_derivation = options_.analytic_derivation;
            options_.analytic_derivation = -2;
            info = stoch_simul(options_.varobs);
            data_info.data=oo_.endo_simul(options_.varobs_id,100+1:end);
            %                         datax=data;
            derivatives_info.no_DLIK=1;
            [fval,DLIK,AHess,cost_flag,ys,trend_coeff,info,M_,options_,bayestopt_,oo_] = dsge_likelihood(params',data_info,options_,M_,estim_params_,bayestopt_,oo_,derivatives_info);
%                 fval = DsgeLikelihood(xparam1,data_info,options_,M_,estim_params_,bayestopt_,oo_);
            options_.analytic_derivation = analytic_derivation;
            AHess=-AHess;
            if min(eig(AHess))<0,
                error('Analytic Hessian is not positive semi-definite!')
            end
%             chol(AHess);
            ide_hess.AHess= AHess;
            deltaM = sqrt(diag(AHess));
            iflag=any((deltaM.*deltaM)==0);
            tildaM = AHess./((deltaM)*(deltaM'));
            if iflag || rank(AHess)>rank(tildaM),
                [ide_hess.cond, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(AHess, 1);
            else
                [ide_hess.cond, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(tildaM, 1);
            end
            indok = find(max(ide_hess.indno,[],1)==0);
            cparam(indok,indok) = inv(AHess(indok,indok));
            normaliz(indok) = sqrt(diag(cparam(indok,indok)))';
            cmm = NaN(size(siJ,1),size(siJ,1));
            ind1=find(ide_hess.ind0);
            cmm = siJ(:,ind1)*((AHess(ind1,ind1))\siJ(:,ind1)');
            temp1=((AHess(ind1,ind1))\siH(:,ind1)');
            diag_chh=sum(siH(:,ind1)'.*temp1)';
%             chh = siH(:,ind1)*((AHess(ind1,ind1))\siH(:,ind1)');
            ind1=ind1(ind1>offset);
            clre = siLRE(:,ind1-offset)*((AHess(ind1,ind1))\siLRE(:,ind1-offset)');
            rhoM=sqrt(1./diag(inv(tildaM(indok,indok))));
%             deltaM = deltaM.*abs(params');
            flag_score=1;
        catch,
            replic = max([replic, length(indJJ)*3]);
            cmm = simulated_moment_uncertainty(indJJ, periods, replic);
%             MIM=siJ(:,indok)'*(cmm\siJ(:,indok));
%           look for independent moments!
            sd=sqrt(diag(cmm));
            cc=cmm./(sd*sd');
            ix=[];
            for jc=1:length(cmm),
                jcheck=find(abs(cc(:,jc))>(1-1.e-6));
                ix=[ix; jcheck(jcheck>jc)];
            end
            iy=find(~ismember([1:length(cmm)],ix));
            indJJ=indJJ(iy);
            GAM=GAM(iy);
            cmm=cmm(iy,iy);
            siJ = (JJ(indJJ,:));
            MIM=siJ'*(cmm\siJ);
            ide_hess.AHess= MIM;
            deltaM = sqrt(diag(MIM));
            iflag=any((deltaM.*deltaM)==0);
            tildaM = MIM./((deltaM)*(deltaM'));
            if iflag || rank(MIM)>rank(tildaM),
                [ide_hess.cond, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(MIM, 1);
            else
                [ide_hess.cond, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(tildaM, 1);
            end
            indok = find(max(ide_hess.indno,[],1)==0);
%             rhoM=sqrt(1-1./diag(inv(tildaM)));
%             rhoM=(1-1./diag(inv(tildaM)));
            ind1=find(ide_hess.ind0);
            temp1=((MIM(ind1,ind1))\siH(:,ind1)');
            diag_chh=sum(siH(:,ind1)'.*temp1)';
%             chh = siH(:,ind1)*((MIM(ind1,ind1))\siH(:,ind1)');
            ind1=ind1(ind1>offset);
            clre = siLRE(:,ind1-offset)*((MIM(ind1,ind1))\siLRE(:,ind1-offset)');
            if ~isempty(indok),
                rhoM(indok)=sqrt(1./diag(inv(tildaM(indok,indok))));
                normaliz(indok) = (sqrt(diag(inv(tildaM(indok,indok))))./deltaM(indok))'; %sqrt(diag(inv(MIM(indok,indok))))';
            end
            %             deltaM = deltaM.*abs(params')
            flag_score=0;
        end
        ide_strength_J(indok) = (1./(normaliz(indok)'./abs(params(indok)')));
        ide_strength_J_prior(indok) = (1./(normaliz(indok)'./normaliz1(indok)'));
        ide_strength_J(params==0)=ide_strength_J_prior(params==0);
        deltaM_prior = deltaM.*abs(normaliz1');
        deltaM = deltaM.*abs(params');
        deltaM(params==0)=deltaM_prior(params==0);
        quant = siJ./repmat(sqrt(diag(cmm)),1,nparam);
        siJnorm = vnorm(quant).*normaliz1;
        %                 siJnorm = vnorm(siJ(inok,:)).*normaliz;
        quant=[];
%         inok = find((abs(TAU)<1.e-8));
%         isok = find((abs(TAU)>=1.e-8));
%         quant(isok,:) = siH(isok,:)./repmat(TAU(isok,1),1,nparam);
%         quant(inok,:) = siH(inok,:)./repmat(mean(abs(TAU)),length(inok),nparam);
%         quant = siH./repmat(sqrt(diag(chh)),1,nparam);
        iy = find(diag_chh);
        indH=indH(iy);
        siH=siH(iy,:);
        if ~isempty(iy),
            quant = siH./repmat(sqrt(diag_chh(iy)),1,nparam);
            siHnorm = vnorm(quant).*normaliz1;
        else
            siHnorm = [];
        end
        %                 siHnorm = vnorm(siH./repmat(TAU,1,nparam)).*normaliz;
        quant=[];
%         inok = find((abs(LRE)<1.e-8));
%         isok = find((abs(LRE)>=1.e-8));
%         quant(isok,:) = siLRE(isok,:)./repmat(LRE(isok,1),1,np);
%         quant(inok,:) = siLRE(inok,:)./repmat(mean(abs(LRE)),length(inok),np);
        diag_clre = diag(clre);
        iy = find(diag_clre);
        indLRE=indLRE(iy);
        siLRE=siLRE(iy,:);
        if ~isempty(iy),
            quant = siLRE./repmat(sqrt(diag_clre(iy)),1,np);
            siLREnorm = vnorm(quant).*normaliz1(offset+1:end);
        else
            siLREnorm=[];
        end
        %                 siLREnorm = vnorm(siLRE./repmat(LRE,1,nparam-offset)).*normaliz(offset+1:end);
        ide_hess.ide_strength_J=ide_strength_J; 
        ide_hess.ide_strength_J_prior=ide_strength_J_prior; 
        ide_hess.deltaM=deltaM; 
        ide_hess.deltaM_prior=deltaM_prior; 
        ide_moments.siJnorm=siJnorm; 
        ide_model.siHnorm=siHnorm; 
        ide_lre.siLREnorm=siLREnorm; 
        ide_hess.flag_score=flag_score; 
    end,
    if normalize_jacobians,
        normH = max(abs(siH)')';
        normH = normH(:,ones(nparam,1));
        normJ = max(abs(siJ)')';
        normJ = normJ(:,ones(nparam,1));
        normLRE = max(abs(siLRE)')';
        normLRE = normLRE(:,ones(size(gp,2),1));
    else
        normH = 1;
        normJ = 1;
        normLRE = 1;
    end
    ide_moments.indJJ=indJJ;
    ide_model.indH=indH;
    ide_lre.indLRE=indLRE;
    ide_moments.siJ=siJ;
    ide_model.siH=siH;
    ide_lre.siLRE=siLRE;
    ide_moments.GAM=GAM;
    ide_model.TAU=TAU;
    ide_lre.LRE=LRE;
%     [ide_checks.idemodel_Mco, ide_checks.idemoments_Mco, ide_checks.idelre_Mco, ...
%         ide_checks.idemodel_Pco, ide_checks.idemoments_Pco, ide_checks.idelre_Pco, ...
%         ide_checks.idemodel_cond, ide_checks.idemoments_cond, ide_checks.idelre_cond, ...
%         ide_checks.idemodel_ee, ide_checks.idemoments_ee, ide_checks.idelre_ee, ...
%         ide_checks.idemodel_ind, ide_checks.idemoments_ind, ...
%         ide_checks.idemodel_indno, ide_checks.idemoments_indno, ...
%         ide_checks.idemodel_ino, ide_checks.idemoments_ino] = ...
%         identification_checks(H(indH,:)./normH(:,ones(nparam,1)),JJ(indJJ,:)./normJ(:,ones(nparam,1)), gp(indLRE,:)./normLRE(:,ones(size(gp,2),1)));
    [ide_moments.cond, ide_moments.ind0, ide_moments.indno, ide_moments.ino, ide_moments.Mco, ide_moments.Pco, ide_moments.jweak, ide_moments.jweak_pair] = ...
        identification_checks(JJ(indJJ,:)./normJ, 0);
    [ide_model.cond, ide_model.ind0, ide_model.indno, ide_model.ino, ide_model.Mco, ide_model.Pco, ide_model.jweak, ide_model.jweak_pair] = ...
        identification_checks(H(indH,:)./normH, 0);
    [ide_lre.cond, ide_lre.ind0, ide_lre.indno, ide_lre.ino, ide_lre.Mco, ide_lre.Pco, ide_lre.jweak, ide_lre.jweak_pair] = ...
        identification_checks(gp(indLRE,:)./normLRE, 0);
    normJ=1;
    [U, S, V]=svd(JJ(indJJ,:)./normJ,0);
    S=diag(S);
    S=[S;zeros(size(JJ,2)-length(indJJ),1)];
    if nparam>8
        ide_moments.S = S([1:4, end-3:end]);
        ide_moments.V = V(:,[1:4, end-3:end]);
    else
        ide_moments.S = S;
        ide_moments.V = V;
    end
    
    indok = find(max(ide_moments.indno,[],1)==0);
    if advanced,
        [ide_moments.pars, ide_moments.cosnJ] = ident_bruteforce(JJ(indJJ,:)./normJ,max_dim_cova_group,options_.TeX,name_tex);
    end
end    
