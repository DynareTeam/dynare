function [H, dA, dOm, Hss, gp, d2A, d2Om, H2ss] = getH(A, B, estim_params_,M_,oo_,options_,kronflag,indx,indexo,iv)
% function [H, dA, dOm, Hss, gp, d2A, d2Om, H2ss] = getH(A, B, estim_params_,M_,oo_,options_,kronflag,indx,indexo,iv)
% computes derivative of reduced form linear model w.r.t. deep params
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
%   iv:             Index of considered variables           
% 
% Outputs:
%   H:              dTAU/dTHETA: Jacobian of TAU, vectorized form of
%                   linearized reduced form state space model, given ys [steady state],
%                   A [transition matrix], B [matrix of shocks], Sigma [covariance of shocks]
%                   TAU = [ys; vec(A); dyn_vech(B*Sigma*B')]. 
%   dA:             [endo_nbr by endo_nbr by (indx+indexo)] Jacobian of transition matrix A            
%   dOm:            [endo_nbr by endo_nbr by (indx+indexo)] Jacobian of Omega = (B*Sigma*B')
%   Hss:            [endo_nbr by (indx)] Jacobian of steady state with respect to estimated
%                   structural parameters only (indx)
%   gp:             Jacobian of linear rational expectation matrices [i.e.
%                   Jacobian of dynamic model] with respect to estimated
%                   structural parameters only (indx)
%   d2A:            Hessian of transition matrix A
%   d2Om:           Hessian of Omega
%   H2s:            Hessian of steady state with respect to estimated
%                   structural parameters only (indx)

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

if nargin<7 || isempty(kronflag)
    kronflag = 0; 
end
if nargin<8 || isempty(indx)
    indx = []; 
end
if nargin<9 || isempty(indexo)
    indexo = []; 
end
if nargin<10 || isempty(iv)
    iv = (1:length(A))'; 
end

[I,J]=find(M_.lead_lag_incidence');
yy0=oo_.dr.ys(I);
param_nbr = length(indx);
tot_param_nbr = param_nbr + length(indexo);
if nargout>5,
    param_nbr_2 = param_nbr*(param_nbr+1)/2;
    tot_param_nbr_2 = tot_param_nbr*(tot_param_nbr+1)/2;
end

m = size(A,1);
m1=length(iv);
n = size(B,2);

if kronflag==-1, % perturbation
    gp=0;
    fun = 'thet2tau';
    params0 = M_.params;
    H = fjaco(fun,[sqrt(diag(M_.Sigma_e(indexo,indexo))); M_.params(indx)], M_, oo_, indx, indexo,0);
    if nargout>1,
        dOm = zeros(m1,m1,tot_param_nbr);
        dA=zeros(m1,m1,tot_param_nbr);
        Hss=H(iv,length(indexo)+1:end);
        da = H(m+1:m+m*m,:);
        dom = H(m+m*m+1:end,:);
        for j=1:tot_param_nbr,
            tmp = dyn_unvech(dom(:,j));
            dOm(:,:,j) = tmp(iv,iv);
            tmp = reshape(da(:,j),m,m);
            dA(:,:,j) = tmp(iv,iv);
        end
        clear da dom tmp
    end
    if nargout>5,
        H2 = hessian_sparse('thet2tau',[sqrt(diag(M_.Sigma_e(indexo,indexo))); M_.params(indx)], ...
            options_.gstep,estim_params_,M_, oo_, indx,indexo,0,[],[],[],iv);
        H2ss = zeros(m1,tot_param_nbr,tot_param_nbr);
        iax=find(triu(rand(tot_param_nbr,tot_param_nbr)));
        H2 = H2(:,iax);
        for j=1:m1,
            H2ss(j,:,:)=dyn_unvech(full(H2(j,:)));
        end
        H2ss=H2ss(:,length(indexo)+1:end,length(indexo)+1:end);
        d2A = sparse(m1*m1,tot_param_nbr_2);
        d2Om = sparse(m1*(m1+1)/2,tot_param_nbr_2);
        d2A(:,:) = H2(m1+1:m1+m1*m1,:);
        d2Om(:,:) = H2(m1+m1*m1+1:end,:);
        clear H2
%         tmp0=zeros(m1,m1);
%         tmp0(iv,iv)=1;
%         iax=find(tmp0);
%         d2A=d2a(iax,:);
%         iax=find(dyn_vech(tmp0));
%         d2Om=d2om(iax,:);

    end
%     assignin('base','M_', M_);
%     assignin('base','oo_', oo_);
    return
end

if kronflag==-2,
    if nargout>5,
        [residual, g1, g2 ] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
            M_.params, oo_.dr.ys, 1);
        g22 = hessian_sparse('thet2tau',[M_.params(indx)],options_.gstep,estim_params_,M_, oo_, indx,[],-1);
        H2ss=full(g22(1:M_.endo_nbr,:));
        H2ss = reshape(H2ss,[M_.endo_nbr param_nbr param_nbr]);
        for j=1:M_.endo_nbr,
            H2ss(j,:,:)=dyn_unvech(dyn_vech(H2ss(j,:,:)));
        end
        g22=g22(M_.endo_nbr+1:end,:);
        inx=find(g22);
        gx22=zeros(length(inx),5);
        for j=1:length(inx),
            [i1, i2] = ind2sub(size(g22),inx(j));
            [ig1, ig2] = ind2sub(size(g1),i1);
            [ip1, ip2] = ind2sub([param_nbr param_nbr],i2);
            gx22(j,:) = [ig1 ig2 ip1 ip2 g22(inx(j))];
        end
        g22 = gx22;
        clear gx22;
    else
        [residual, g1 ] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
            M_.params, oo_.dr.ys, 1);        
    end
    gp = fjaco('thet2tau',[M_.params(indx)],estim_params_,M_, oo_, indx,[],-1);
    Hss=gp(1:M_.endo_nbr,:);
    gp=gp(M_.endo_nbr+1:end,:);
    gp = reshape(gp,[size(g1) param_nbr]);
else

% yy0=[];
% for j=1:size(M_.lead_lag_incidence,1);
%     yy0 = [ yy0; oo_.dr.ys(find(M_.lead_lag_incidence(j,:)))];
% end
dyssdtheta=zeros(length(oo_.dr.ys),M_.param_nbr);
d2yssdtheta=zeros(length(oo_.dr.ys),M_.param_nbr,M_.param_nbr);
[residual, gg1] = feval([M_.fname,'_static'],oo_.dr.ys, oo_.exo_steady_state', M_.params);
df = feval([M_.fname,'_static_params_derivs'],oo_.dr.ys, repmat(oo_.exo_steady_state',[M_.maximum_exo_lag+M_.maximum_exo_lead+1]), ...
    M_.params);
dyssdtheta = -gg1\df;
if nargout>5,
    [residual, gg1, gg2] = feval([M_.fname,'_static'],oo_.dr.ys, oo_.exo_steady_state', M_.params);
    [residual, g1, g2, g3] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1);
    [nr, nc]=size(gg2);

    [df, gpx, d2f] = feval([M_.fname,'_static_params_derivs'],oo_.dr.ys, oo_.exo_steady_state', ...
        M_.params);%, oo_.dr.ys, 1, dyssdtheta*0, d2yssdtheta);
    d2f = get_all_resid_2nd_derivs(d2f,length(oo_.dr.ys),M_.param_nbr);

    if isempty(find(gg2)),
        for j=1:M_.param_nbr,
        d2yssdtheta(:,:,j) = -gg1\d2f(:,:,j);
        end
    else
        gam = d2f*0;
        for j=1:nr,
            tmp1 = (squeeze(gpx(j,:,:))'*dyssdtheta);
            gam(j,:,:)=transpose(reshape(gg2(j,:),[nr nr])*dyssdtheta)*dyssdtheta ...
                + tmp1 + tmp1';
        end
        for j=1:M_.param_nbr,
        d2yssdtheta(:,:,j) = -gg1\(d2f(:,:,j)+gam(:,:,j));
%         d2yssdtheta(:,:,j) = -gg1\(d2f(:,:,j)+gam(:,:,j)+ squeeze(gpx(:,:,j))*dyssdtheta);
        end
        clear tmp1 gpx gam,
    end
end

if any(any(isnan(dyssdtheta))),    
    [U,T] = schur(gg1);
    qz_criterium=options_.qz_criterium;
    e1 = abs(ordeig(T)) < qz_criterium-1;
    k = sum(e1);       % Number non stationary variables.
%     n = length(e1)-k;  % Number of stationary variables.
    [U,T] = ordschur(U,T,e1);
    T = T(k+1:end,k+1:end);
    dyssdtheta = -U(:,k+1:end)*(T\U(:,k+1:end)')*df;
    if nargout>5,
        for j=1:length(indx),
            d2yssdtheta(:,:,j) = -U(:,k+1:end)*(T\U(:,k+1:end)')*d2f(:,:,j);
        end
    end
end
if nargout>5,
    [df, gp, d2f, gpp, hp] = feval([M_.fname,'_params_derivs'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1, dyssdtheta, d2yssdtheta);
    H2ss = d2yssdtheta(oo_.dr.order_var,indx,indx);
%     nelem=size(g1,2);
%     g22 = get_all_2nd_derivs(gpp,m,nelem,M_.param_nbr);
%     g22 = g22(:,:,indx,indx);
else
    [df, gp] = feval([M_.fname,'_params_derivs'],yy0, repmat(oo_.exo_steady_state',[M_.maximum_exo_lag+M_.maximum_exo_lead+1,1]), ...
        M_.params, oo_.dr.ys, 1, dyssdtheta,d2yssdtheta);
    [residual, g1, g2 ] = feval([M_.fname,'_dynamic'],yy0, repmat(oo_.exo_steady_state',[M_.maximum_exo_lag+M_.maximum_exo_lead+1,1]), ...
        M_.params, oo_.dr.ys, 1);
end

[nr, nc]=size(g2);
nc = sqrt(nc);
Hss = dyssdtheta(oo_.dr.order_var,indx);
dyssdtheta = dyssdtheta(I,:);
ns = max(max(M_.lead_lag_incidence)); % retrieve the number of states excluding columns for shocks
gp2 = gp*0;
for j=1:nr,
    [II JJ]=ind2sub([nc nc],find(g2(j,:)));
    for i=1:nc,
        is = find(II==i);
        is = is(find(JJ(is)<=ns));
        if ~isempty(is),
            g20=full(g2(j,find(g2(j,:))));
            gp2(j,i,:)=g20(is)*dyssdtheta(JJ(is),:);
        end
    end
end

gp = gp+gp2;
gp = gp(:,:,indx);

if nargout>5,
%     h22 = get_all_hess_derivs(hp,nr,nc,M_.param_nbr);
    g22 = gpp;
    gp22 = sparse(nr*nc,param_nbr*param_nbr);
    tmp1 = reshape(g3,[nr*nc*nc nc]);
    tmp2=sparse(size(tmp1,1),M_.param_nbr);
%     tmp2=tmp1*[dyssdtheta; zeros(nc-ns,M_.param_nbr)];
    tmpa=[dyssdtheta; zeros(nc-ns,M_.param_nbr)];
    tmpa=sparse(tmpa);
    for j=1:M_.param_nbr,
        tmp2(:,j)=tmp1*tmpa(:,j);
    end
%     tmp2=sparse(tmp2);
%     [i1 i2]=ind2sub([nc M_.param_nbr],[1:nc*M_.param_nbr]');

    for j=1:nr,
        tmp0=reshape(g2(j,:),[nc nc]);   
        tmp0 = tmp0(:,1:ns)*reshape(d2yssdtheta(I,:,:),[ns,M_.param_nbr*M_.param_nbr]);
        for i=1:nc,
            indo = sub2ind([nr nc nc], ones(nc,1)*j ,ones(nc,1)*i, (1:nc)');
            tmpx = (tmp2(indo,:))'*[dyssdtheta; zeros(nc-ns,M_.param_nbr)];
%             gp22(j,i,:,:)=squeeze(tmp1(j,i,:,:))'*[dyssdtheta; zeros(nc-ns,M_.param_nbr)];
            tmpu = (get_hess_deriv(hp,j,i,nc,M_.param_nbr))'*[dyssdtheta; zeros(nc-ns,M_.param_nbr)];
            tmpy = tmpx+tmpu+tmpu'+reshape(tmp0(i,:,:),[M_.param_nbr M_.param_nbr]);
            tmpy = tmpy + get_2nd_deriv_mat(gpp,j,i,M_.param_nbr);
            tmpy = tmpy(indx,indx);
            if any(any(tmpy)),
                ina = find(triu(tmpy));
                gp22(sub2ind([nr nc],j,i),ina)=transpose(tmpy(ina));
%             gp22(j,i,:,:)= reshape(tmpy,[1 1 M_.param_nbr M_.param_nbr]);

            end
        end
%             gp22(j,:,:,:)=gp22(j,:,:,:)+reshape(tmp0(:,1:ns)*d2yssdtheta(I,:,:),[1 nc M_.param_nbr M_.param_nbr]);
    end

%     g22 = g22+gp22;
%     g22 = g22(:,:,indx,indx);
    clear tmp0 tmp1 tmp2 tmpu tmpx tmpy,
        inx=find(gp22);
        gx22=zeros(length(inx),5);
        for j=1:length(inx),
            [i1, i2] = ind2sub(size(gp22),inx(j));
            [ig1, ig2] = ind2sub(size(g1),i1);
            [ip1, ip2] = ind2sub([param_nbr param_nbr],i2);
            gx22(j,:) = [ig1 ig2 ip1 ip2 gp22(inx(j))];
        end
        g22 = gx22;
        clear gx22 gp22;

end
end



klen = M_.maximum_endo_lag + M_.maximum_endo_lead + 1;
k11 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_endo_lag+1),:);
a = g1(:,nonzeros(k11'));
da = gp(:,nonzeros(k11'),:);
if nargout > 5,
    indind = ismember(g22(:,2),nonzeros(k11'));
    tmp = g22(indind,:);
    d2a=tmp;
    for j=1:size(tmp,1),
        inxinx = find(nonzeros(k11')==tmp(j,2));
        d2a(j,2) = inxinx;
    end
%     d2a = g22(:,nonzeros(k11'),:,:);
end
kstate = oo_.dr.kstate;

GAM1 = zeros(M_.endo_nbr,M_.endo_nbr);
Dg1 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr);

k1 = find(kstate(:,2) == M_.maximum_endo_lag+2 & kstate(:,3));
GAM1(:, kstate(k1,1)) = -a(:,kstate(k1,3));
Dg1(:, kstate(k1,1), :) = -da(:,kstate(k1,3),:);
if nargout > 5,
    indind = ismember(d2a(:,2),kstate(k1,3));
    tmp = d2a(indind,:);
    tmp(:,end)=-tmp(:,end);
    D2g1 = tmp;
    for j=1:size(tmp,1),
        inxinx = (kstate(k1,3)==tmp(j,2));
        D2g1(j,2) = kstate(k1(inxinx),1);
    end
end

[junk,cols_b,cols_j] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+1, ...
                                                  oo_.dr.order_var));
GAM0 = zeros(M_.endo_nbr,M_.endo_nbr);
Dg0 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr);
GAM0(:,cols_b) = g1(:,cols_j);
Dg0(:,cols_b,:) = gp(:,cols_j,:);
if nargout > 5,
    indind = ismember(g22(:,2),cols_j);
    tmp = g22(indind,:);
    D2g0=tmp;
    for j=1:size(tmp,1),
        inxinx = (cols_j==tmp(j,2));
        D2g0(j,2) = cols_b(inxinx);
    end
end


k2 = find(kstate(:,2) == M_.maximum_endo_lag+1 & kstate(:,4));
GAM2 = zeros(M_.endo_nbr,M_.endo_nbr);
Dg2 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr);
GAM2(:, kstate(k2,1)) = -a(:,kstate(k2,4));
Dg2(:, kstate(k2,1), :) = -da(:,kstate(k2,4),:);
if nargout > 5,
    indind = ismember(d2a(:,2),kstate(k2,4));
    tmp = d2a(indind,:);
    tmp(:,end)=-tmp(:,end);
    D2g2 = tmp;
    for j=1:size(tmp,1),
        inxinx = (kstate(k2,4)==tmp(j,2));
        D2g2(j,2) = kstate(k2(inxinx),1);
    end
end

GAM3 = -g1(:,length(yy0)+1:end);
Dg3 = -gp(:,length(yy0)+1:end,:);
if nargout>5,
    cols_ex = [length(yy0)+1:size(g1,2)];
    indind = ismember(g22(:,2),cols_ex);
    tmp = g22(indind,:);
    tmp(:,end)=-tmp(:,end);
    D2g3=tmp;
    for j=1:size(tmp,1),
        inxinx = find(cols_ex==tmp(j,2));
        D2g3(j,2) = inxinx;
    end
    clear g22 d2a tmp
end

clear g1 g2 g3 df d2f gpp hp residual gg1 gg2 gp2 dyssdtheta d2yssdtheta

if kronflag==1, % kronecker products
    Dg0=reshape(Dg0,m^2,param_nbr);
    Dg1=reshape(Dg1,m^2,param_nbr);
    Dg2=reshape(Dg2,m^2,param_nbr);
    for j=1:param_nbr,
        Dg3(:,:,j)=Dg3(:,:,j)*M_.Sigma_e;
    end
    Dg3=reshape(Dg3,m*n,param_nbr);
    Om = B*M_.Sigma_e*B';
    Im = eye(m);
    Dm = duplication(m);
    DmPl = inv(Dm'*Dm)*Dm';
    Kmm = commutation(m,m);
    Kmn = commutation(m,n);


    Da = [eye(m^2),zeros(m^2,m*(m+1)/2)];
    Dom = [zeros(m*(m+1)/2,m^2),eye(m*(m+1)/2)];


    Df1Dtau = ( kron(Im,GAM0) - kron(A',GAM1) - kron(Im,GAM1*A) )*Da;

    Df1Dthet = kron(A',Im)*Dg0 - kron( (A')^2,Im)*Dg1 - Dg2;

    Df2Dtau = DmPl*( kron(GAM0,GAM0) - kron(GAM0,GAM1*A) - kron(GAM1*A,GAM0) + kron(GAM1*A,GAM1*A) )*Dm*Dom - ...
              DmPl*( kron(GAM0*Om,GAM1) + kron(GAM1,GAM0*Om)*Kmm - kron(GAM1*A*Om,GAM1) - kron(GAM1,GAM1*A*Om)*Kmm )*Da;


    Df2Dthet = DmPl*( kron(GAM0*Om,Im) + kron(Im,GAM0*Om)*Kmm - kron(Im,GAM1*A*Om)*Kmm - kron(GAM1*A*Om,Im) )*Dg0 - ...
        DmPl*( kron(GAM0*Om*A',Im) + kron(Im,GAM0*Om*A')*Kmm - kron(Im,GAM1*A*Om*A')*Kmm - kron(GAM1*A*Om*A',Im) )*Dg1 -...
        DmPl*( kron(GAM3,Im) + kron(Im,GAM3)*Kmn )*Dg3;


    DfDtau  = [Df1Dtau;Df2Dtau];
    DfDthet = [Df1Dthet;Df2Dthet];

    H = -DfDtau\DfDthet;
    x = reshape(H(1:m*m,:),m,m,param_nbr);
    y = reshape(Dm*H(m*m+1:end,:),m,m,param_nbr);
    dA = x;
    dOm = y;
    % convert to dyn_vech
    tmpH = Dm*H(m*m+1:end,:);
    Index = find(triu(ones(m)));
    H(m*m+1:end,:) = tmpH(Index,:);

    Hx = [];
    if ~isempty(indexo),
        dSig = zeros(M_.exo_nbr,M_.exo_nbr);
        dOm = cat(3,zeros(size(dOm,1),size(dOm,1),length(indexo)),dOm);
        for j=1:length(indexo)
            dSig(indexo(j),indexo(j)) = 2*sqrt(M_.Sigma_e(indexo(j),indexo(j)));
            y = B*dSig*B';
            y = y(nauxe+1:end,nauxe+1:end);
            Hx(:,j) = [zeros((m-nauxe)^2,1); dyn_vech(y)];
            if nargout>1,
                dOm(:,:,j) = y;
            end
            dSig(indexo(j),indexo(j)) = 0;
        end
    end
    H = [ [zeros(M_.endo_nbr,length(indexo)) Hss]; [Hx H]];

else % generalized sylvester equation

    % solves a*x+b*x*c=d
    a = (GAM0-GAM1*A);
    inva = inv(a);
    b = -GAM1;
    c = A;
    elem = zeros(m,m,param_nbr);
    d = elem;
    for j=1:param_nbr,
        elem(:,:,j) = (Dg0(:,:,j)-Dg1(:,:,j)*A);
        d(:,:,j) = Dg2(:,:,j)-elem(:,:,j)*A;
    end
    xx=sylvester3(a,b,c,d);
    flag=1;
    icount=0;
    while flag && icount<4,
        [xx, flag]=sylvester3a(xx,a,b,c,d);
        icount=icount+1;
    end
    H=zeros(m1*m1+m1*(m1+1)/2,param_nbr+length(indexo));
    if nargout>1,
        dOm = zeros(m1,m1,param_nbr+length(indexo));
        dA=zeros(m1,m1,param_nbr+length(indexo));
        dB=zeros(m,n,param_nbr);
    end
    if ~isempty(indexo),
        dSig = zeros(M_.exo_nbr,M_.exo_nbr,length(indexo));
        for j=1:length(indexo)
            dSig(indexo(j),indexo(j),j) = 2*sqrt(M_.Sigma_e(indexo(j),indexo(j)));
            y = B*dSig(:,:,j)*B';
%             y = y(nauxe+1:end,nauxe+1:end);
%             H(:,j) = [zeros((m-nauxe)^2,1); dyn_vech(y)];
            H(:,j) = [zeros(m1^2,1); dyn_vech(y(iv,iv))];
            if nargout>1,
                dOm(:,:,j) = y(iv,iv);
            end
%             dSig(indexo(j),indexo(j)) = 0;
        end
    end
    for j=1:param_nbr,
        x = xx(:,:,j);
        y = inva * (Dg3(:,:,j)-(elem(:,:,j)-GAM1*x)*B);
        if nargout>1,
            dB(:,:,j) = y;
        end
        y = y*M_.Sigma_e*B'+B*M_.Sigma_e*y';
%         x = x(nauxe+1:end,nauxe+1:end);
%         y = y(nauxe+1:end,nauxe+1:end);
        if nargout>1,
            dA(:,:,j+length(indexo)) = x(iv,iv);
            dOm(:,:,j+length(indexo)) = y(iv,iv);
        end
        H(:,j+length(indexo)) = [vec(x(iv,iv)); dyn_vech(y(iv,iv))];
    end
    %     for j=1:param_nbr,
    %         disp(['Derivatives w.r.t. ',M_.param_names(indx(j),:),', ',int2str(j),'/',int2str(param_nbr)])
    %         elem = (Dg0(:,:,j)-Dg1(:,:,j)*A);
    %         d = Dg2(:,:,j)-elem*A;
    %         x=sylvester3(a,b,c,d);
    % %         x=sylvester3a(x,a,b,c,d);
    %         y = inva * (Dg3(:,:,j)-(elem-GAM1*x)*B);
    %         y = y*B'+B*y';
    %         x = x(nauxe+1:end,nauxe+1:end);
    %         y = y(nauxe+1:end,nauxe+1:end);
    %         H(:,j) = [x(:); dyn_vech(y)];
    %     end
    Hss = Hss(iv,:);
    H = [[zeros(m1,length(indexo)) Hss]; H];

end

if nargout > 5,
    H2ss = H2ss(iv,:,:);
    d = zeros(m,m,floor(sqrt(param_nbr_2)));
    %     d2A = zeros(m,m,tot_param_nbr,tot_param_nbr);
    %     d2Om = zeros(m,m,tot_param_nbr,tot_param_nbr);
    %     d2B = zeros(m,n,tot_param_nbr,tot_param_nbr);
    %     cc=triu(ones(param_nbr,param_nbr));
    %     [i2,j2]=find(cc);
    %     cc = blkdiag( zeros(length(indexo),length(indexo)), cc);
    %     [ipar2]=find(cc);
    %     ctot=triu(ones(tot_param_nbr,tot_param_nbr));
    %     ctot(1:length(indexo),1:length(indexo))=eye(length(indexo));
    %     [itot2, jtot2]=find(ctot);
    jcount=0;
    cumjcount=0;
    jinx = [];
    x2x=sparse(m*m,param_nbr_2);
%     x2x=[];
    for i=1:param_nbr,
        for j=1:i,
            elem1 = (get_2nd_deriv(D2g0,m,m,j,i)-get_2nd_deriv(D2g1,m,m,j,i)*A);
            elem1 = get_2nd_deriv(D2g2,m,m,j,i)-elem1*A;
            elemj0 = Dg0(:,:,j)-Dg1(:,:,j)*A;
            elemi0 = Dg0(:,:,i)-Dg1(:,:,i)*A;
            elem2 = -elemj0*xx(:,:,i)-elemi0*xx(:,:,j);
            elem2 = elem2 + ( Dg1(:,:,j)*xx(:,:,i) + Dg1(:,:,i)*xx(:,:,j) )*A;
            elem2 = elem2 + GAM1*( xx(:,:,i)*xx(:,:,j) + xx(:,:,j)*xx(:,:,i));
            jcount=jcount+1;
            jinx = [jinx; [j i]];
            d(:,:,jcount) = elem1+elem2;
            if jcount==floor(sqrt(param_nbr_2)) || (j*i)==param_nbr^2,
                if (j*i)==param_nbr^2,
                    d = d(:,:,1:jcount);
                end
%                 d(find(abs(d)<1.e-12))=0;
                xx2=sylvester3(a,b,c,d);
                flag=1;
                icount=0;
                while flag && icount<4,
                    [xx2, flag]=sylvester3a(xx2,a,b,c,d);
                    icount = icount + 1;
                end
%                 inx = find(abs(xx2)>1.e-12);
%                 xx2(find(abs(xx2)<1.e-12))=0;
                x2x(:,cumjcount+1:cumjcount+jcount)=reshape(xx2,[m*m jcount]);
                cumjcount=cumjcount+jcount;
%                 [i1 i2 i3]=ind2sub(size(xx2),inx);
%                 x2x = [x2x; [i1 i2 jinx(i3,:) xx2(inx)]];
                jcount = 0;
                jinx = [];
            end
        end
    end
    clear d xx2;
    jcount = 0;
    icount = 0;
    cumjcount = 0;
    MAX_DIM_MAT = 100000000;
    ncol = max(1,floor(MAX_DIM_MAT/(8*m1*(m1+1)/2)));
    ncol = min(ncol, tot_param_nbr_2);
    d2A = sparse(m1*m1,tot_param_nbr_2);
    d2Om = sparse(m1*(m1+1)/2,tot_param_nbr_2);
    d2A_tmp = zeros(m1*m1,ncol);
    d2Om_tmp = zeros(m1*(m1+1)/2,ncol);
    tmpDir = CheckPath('tmp_derivs',M_.dname);
    offset = length(indexo);
    %     d2B = zeros(m,n,tot_param_nbr,tot_param_nbr);
    d2Sig = zeros(M_.exo_nbr,M_.exo_nbr,length(indexo));
    for j=1:tot_param_nbr,
        for i=1:j,
            jcount=jcount+1;
            if j<=offset,
                if i==j,
                    d2Sig(indexo(j),indexo(j),j) = 2;
                    y = B*d2Sig(:,:,j)*B';
%                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y(iv,iv));
                end
            else
                jind = j-offset;
                iind = i-offset;
                if i<=offset,
                    y = dB(:,:,jind)*dSig(:,:,i)*B'+B*dSig(:,:,i)*dB(:,:,jind)';
%                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y(iv,iv));
                else
                    icount=icount+1;
                    x = reshape(x2x(:,icount),[m m]);
%                     x = get_2nd_deriv(x2x,m,m,iind,jind);%xx2(:,:,jcount);
                    elem1 = (get_2nd_deriv(D2g0,m,m,iind,jind)-get_2nd_deriv(D2g1,m,m,iind,jind)*A);
                    elem1 = elem1 -( Dg1(:,:,jind)*xx(:,:,iind) + Dg1(:,:,iind)*xx(:,:,jind) );
                    elemj0 = Dg0(:,:,jind)-Dg1(:,:,jind)*A-GAM1*xx(:,:,jind);
                    elemi0 = Dg0(:,:,iind)-Dg1(:,:,iind)*A-GAM1*xx(:,:,iind);
                    elem0 = elemj0*dB(:,:,iind)+elemi0*dB(:,:,jind);
                    y = inva * (get_2nd_deriv(D2g3,m,n,iind,jind)-elem0-(elem1-GAM1*x)*B);
                    %         d2B(:,:,j+length(indexo),i+length(indexo)) = y;
                    %         d2B(:,:,i+length(indexo),j+length(indexo)) = y;
                    y = y*M_.Sigma_e*B'+B*M_.Sigma_e*y'+ ...
                        dB(:,:,jind)*M_.Sigma_e*dB(:,:,iind)'+dB(:,:,iind)*M_.Sigma_e*dB(:,:,jind)';
%                     x(abs(x)<1.e-8)=0;
                    d2A_tmp(:,jcount) = vec(x(iv,iv));
%                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y(iv,iv));
                end
            end
            if jcount==ncol || i*j==tot_param_nbr^2,
                d2A(:,cumjcount+1:cumjcount+jcount) = d2A_tmp(:,1:jcount);
                %         d2A(:,:,j+length(indexo),i+length(indexo)) = x;
                %         d2A(:,:,i+length(indexo),j+length(indexo)) = x;
                d2Om(:,cumjcount+1:cumjcount+jcount) = d2Om_tmp(:,1:jcount);
                %         d2Om(:,:,j+length(indexo),i+length(indexo)) = y;
                %         d2Om(:,:,i+length(indexo),j+length(indexo)) = y;
                save([tmpDir filesep 'd2A_' int2str(cumjcount+1) '_' int2str(cumjcount+jcount) '.mat'],'d2A')
                save([tmpDir filesep 'd2Om_' int2str(cumjcount+1) '_'  int2str(cumjcount+jcount) '.mat'],'d2Om')
                cumjcount = cumjcount+jcount;
                jcount=0;
                %         d2A = sparse(m1*m1,tot_param_nbr*(tot_param_nbr+1)/2);
                %         d2Om = sparse(m1*(m1+1)/2,tot_param_nbr*(tot_param_nbr+1)/2);
                d2A_tmp = zeros(m1*m1,ncol);
                d2Om_tmp = zeros(m1*(m1+1)/2,ncol);
                
            end
        end
    end
end

return

function g22 = get_2nd_deriv(gpp,m,n,i,j),

g22=zeros(m,n);
is=find(gpp(:,3)==i);
is=is(find(gpp(is,4)==j));

if ~isempty(is),
    g22(sub2ind([m,n],gpp(is,1),gpp(is,2)))=gpp(is,5)';
end
return

function g22 = get_2nd_deriv_mat(gpp,i,j,n),

g22=zeros(n,n);
is=find(gpp(:,1)==i);
is=is(find(gpp(is,2)==j));

if ~isempty(is),
    g22(sub2ind([n,n],gpp(is,3),gpp(is,4)))=gpp(is,5)';
    g22(sub2ind([n,n],gpp(is,4),gpp(is,3)))=gpp(is,5)';
end
return

function g22 = get_all_2nd_derivs(gpp,m,n,npar,fsparse),

if nargin==4 || isempty(fsparse),
    fsparse=0;
end
if fsparse,
    g22=sparse(m*n,npar*npar);
else
g22=zeros(m,n,npar,npar);
end
% c=ones(npar,npar);
% c=triu(c);
% ic=find(c);

for is=1:length(gpp),
%     d=zeros(npar,npar);
%     d(gpp(is,3),gpp(is,4))=1;
%     indx = find(ic==find(d));
    if fsparse,
        g22(sub2ind([m,n],gpp(is,1),gpp(is,2)),sub2ind([npar,npar],gpp(is,3),gpp(is,4)))=gpp(is,5);
        g22(sub2ind([m,n],gpp(is,1),gpp(is,2)),sub2ind([npar,npar],gpp(is,4),gpp(is,3)))=gpp(is,5);
    else
    g22(gpp(is,1),gpp(is,2),gpp(is,3),gpp(is,4))=gpp(is,5);
    g22(gpp(is,1),gpp(is,2),gpp(is,4),gpp(is,3))=gpp(is,5);
end
end

return

function r22 = get_all_resid_2nd_derivs(rpp,m,npar),

r22=zeros(m,npar,npar);
% c=ones(npar,npar);
% c=triu(c);
% ic=find(c);

for is=1:length(rpp),
%     d=zeros(npar,npar);
%     d(rpp(is,2),rpp(is,3))=1;
%     indx = find(ic==find(d));
    r22(rpp(is,1),rpp(is,2),rpp(is,3))=rpp(is,4);
    r22(rpp(is,1),rpp(is,3),rpp(is,2))=rpp(is,4);
end

return

function h2 = get_all_hess_derivs(hp,r,m,npar),

h2=zeros(r,m,m,npar);

for is=1:length(hp),
    h2(hp(is,1),hp(is,2),hp(is,3),hp(is,4))=hp(is,5);
    h2(hp(is,1),hp(is,3),hp(is,2),hp(is,4))=hp(is,5);
end

return

function h2 = get_hess_deriv(hp,i,j,m,npar),

h2=zeros(m,npar);
is1=find(hp(:,1)==i);
is=is1(find(hp(is1,2)==j));

if ~isempty(is),
    h2(sub2ind([m,npar],hp(is,3),hp(is,4)))=hp(is,5)';
end

is=is1(find(hp(is1,3)==j));

if ~isempty(is),
    h2(sub2ind([m,npar],hp(is,2),hp(is,4)))=hp(is,5)';
end

return
