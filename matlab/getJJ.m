function [JJ, H, gam, gp, dA, dOm, dYss] = getJJ(A, B, estim_params_, M_,oo_,options_,kronflag,indx,indexo,mf,nlags,useautocorr)
% function [JJ, H, gam, gp, dA, dOm, dYss] = getJJ(A, B, estim_params_, M_,oo_,options_,kronflag,indx,indexo,mf,nlags,useautocorr)
% computes derivatives of 1st and 2nd order moments of observables with
% respect to estimated parameters
% 
% Inputs:
%   A:              Transition matrix of lagged states from Kalman filter
%   B:              Matrix in state transition equation mapping shocks today to
%                   states today
%   M_:             structure storing the model information
%   oo_:            structure storing the results
%   options_:       structure storing the options
%   kronflag:       Indicator whether to rely on Kronecker products (1) or
%                   not (-1 or -2)
%   indx:           Index of estimated parameters in M_.params
%   indexo:         Index of estimated standard deviations in M_.exo_names
%   mf:             Index of observed variables
%   nlags:          Number of lags to consider for covariances and
%                   correlations
%   useautocorr:    Indicator on whether to use correlations (1) instead of
%                   covariances (0)
% 
% Outputs:
%   JJ:             Jacobian of 1st and 2nd order moments of observables, i.e.  dgam/dTHETA 
%                   (see below for definition of gam)
%   H:              dTAU/dTHETA: Jacobian of TAU, vectorized form of
%                   linearized reduced form state space model, given ys [steady state],
%                   A [transition matrix], B [matrix of shocks], Sigma [covariance of shocks]
%                   TAU = [ys; vec(A); dyn_vech(B*Sigma*B')].
%   gam:            vector of theoretical moments of observed variables mf [JJ is the Jacobian of gam]. 
%                   gam = [ys(mf); dyn_vech(GAM{1}); vec(GAM{j+1})]; for j=1:ar and where 
%                   GAM is the first output of th_autocovariances
%   gp:             Jacobian of linear rational expectation matrices [i.e.
%                   Jacobian of dynamic model] with respect to estimated
%                   structural parameters only (indx)
%   dA:             [endo_nbr by endo_nbr by (indx+indexo)] Jacobian of transition matrix A           
%   dOm:            Jacobian of Omega = (B*Sigma*B')
%   dYss            Jacobian of steady state with respect to estimated structural parameters only (indx)

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

if nargin<8 || isempty(indx)
%     indx = [1:M_.param_nbr];
end,
if nargin<9 || isempty(indexo)
    indexo = [];
end,
if nargin<11 || isempty(nlags)
    nlags=3; 
end
if nargin<12 || isempty(useautocorr)
    useautocorr=0; 
end

%   if useautocorr,
warning('off','MATLAB:divideByZero')
%   end
if kronflag == -1,
    fun = 'thet2tau';
    params0 = M_.params;
    para0 = get_all_parameters(estim_params_, M_);
    JJ = fjaco(fun,para0,estim_params_,M_, oo_, indx,indexo,1,mf,nlags,useautocorr);
    M_.params = params0;
    params0 = M_.params;
    H = fjaco(fun,para0,estim_params_,M_, oo_, indx,indexo,0,mf,nlags,useautocorr);
    M_.params = params0;
    params0 = M_.params;
    gp = fjaco(fun,para0,estim_params_,M_, oo_, indx,indexo,-1);
    M_.params = params0;
    offset = length(para0)-length(indx);
    gp = gp(:,offset+1:end);
    dYss = H(1:M_.endo_nbr,offset+1:end);
    dA = reshape(H(M_.orig_endo_nbr+[1:numel(A)],:),[size(A),size(H,2)]);
    dOm = dA*0;
    for j=1:size(H,2),
        dOm(:,:,j) = dyn_unvech(H(M_.endo_nbr+numel(A)+1:end,j));
    end
    assignin('base','M_', M_);
    assignin('base','oo_', oo_);
else
    [H, dA, dOm, dYss, gp] = getH(A, B, estim_params_,M_,oo_,options_,kronflag,indx,indexo);
    gp = reshape(gp,size(gp,1)*size(gp,2),size(gp,3));
    gp = [dYss; gp];
    %   if isempty(H),
    %     JJ = [];
    %     GAM = [];
    %     return
    %   end
    m = length(A);
    GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,1,options_.debug);
    k = find(abs(GAM) < 1e-12);
    GAM(k) = 0;
    %   if useautocorr,
    sdy = sqrt(diag(GAM));
    sy = sdy*sdy';
    %   end
    
    %   BB = dOm*0;
    %   for j=1:length(indx),
    %     BB(:,:,j)= dA(:,:,j)*GAM*A'+A*GAM*dA(:,:,j)'+dOm(:,:,j);
    %   end
    %   XX =  lyapunov_symm_mr(A,BB,options_.qz_criterium,options_.lyapunov_complex_threshold,0);
    for j=1:length(indexo),
        dum =  lyapunov_symm(A,dOm(:,:,j),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,2,options_.debug);
        %     dum =  XX(:,:,j);
        k = find(abs(dum) < 1e-12);
        dum(k) = 0;
        if useautocorr
            dsy = 1/2./sdy.*diag(dum);
            dsy = dsy*sdy'+sdy*dsy';
            dum1=dum;
            dum1 = (dum1.*sy-dsy.*GAM)./(sy.*sy);
            dum1 = dum1-diag(diag(dum1))+diag(diag(dum));
            dumm = dyn_vech(dum1(mf,mf));
        else
            dumm = dyn_vech(dum(mf,mf));
        end
        for i=1:nlags,
            dum1 = A^i*dum;
            if useautocorr
                dum1 = (dum1.*sy-dsy.*(A^i*GAM))./(sy.*sy);
            end
            dumm = [dumm; vec(dum1(mf,mf))];
        end
        JJ(:,j) = dumm;
    end
    nexo = length(indexo);
    for j=1:length(indx),
        dum =  lyapunov_symm(A,dA(:,:,j+nexo)*GAM*A'+A*GAM*dA(:,:,j+nexo)'+dOm(:,:,j+nexo),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,2,options_.debug);
        %     dum =  XX(:,:,j);
        k = find(abs(dum) < 1e-12);
        dum(k) = 0;
        if useautocorr
            dsy = 1/2./sdy.*diag(dum);
            dsy = dsy*sdy'+sdy*dsy';
            dum1=dum;
            dum1 = (dum1.*sy-dsy.*GAM)./(sy.*sy);
            dum1 = dum1-diag(diag(dum1))+diag(diag(dum));
            dumm = dyn_vech(dum1(mf,mf));
        else
            dumm = dyn_vech(dum(mf,mf));
        end
        for i=1:nlags,
            dum1 = A^i*dum;
            for ii=1:i,
                dum1 = dum1 + A^(ii-1)*dA(:,:,j+nexo)*A^(i-ii)*GAM;
            end
            if useautocorr
                dum1 = (dum1.*sy-dsy.*(A^i*GAM))./(sy.*sy);
            end
            dumm = [dumm; vec(dum1(mf,mf))];
        end
        JJ(:,j+nexo) = dumm;
    end
    
    JJ = [ [zeros(length(mf),nexo) dYss(mf,:)]; JJ];
    
end
if nargout >2,
    %     sy=sy(mf,mf);
    options_.ar=nlags;
    nodecomposition = 1;
    [GAM,stationary_vars] = th_autocovariances(oo_.dr,oo_.dr.order_var(mf),M_,options_,nodecomposition);
    sy=sqrt(diag(GAM{1}));
    sy=sy*sy';
    if useautocorr,
        sy=sy-diag(diag(sy))+eye(length(mf));
        GAM{1}=GAM{1}./sy;
    else
        for j=1:nlags,
            GAM{j+1}=GAM{j+1}.*sy;
        end
    end
    gam = dyn_vech(GAM{1});
    for j=1:nlags,
        gam = [gam; vec(GAM{j+1})];
    end
end
gam = [oo_.dr.ys(oo_.dr.order_var(mf)); gam];

%   if useautocorr,
warning('on','MATLAB:divideByZero')
%   end
