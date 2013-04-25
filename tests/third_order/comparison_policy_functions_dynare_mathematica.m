%read in the FV et al. policy functions derived from Mathematica
load policyfunctions

order=options_.order;
dr=oo_.dr;
nx =size(dr.ghx,2);
nu =size(dr.ghu,2);
k = find(dr.kstate(:,2) <= M_.maximum_lag+1);
klag = dr.kstate(k,[1 2]);
state_vars=dr.order_var(klag(:,1));
%order of endogenous variables in FV et al. paper: c, invest, y, h, r, dp, kp, lambda, pi
varlist_FV={'C';'I';'Y';'H';'r';'D';'K';'lambda';'phi'};
for ii=1: length(varlist_FV)
    FV_endo_order(ii,1)=strmatch(varlist_FV(ii),M_.endo_names(dr.order_var,:),'exact');
end

% order in states of FV et al. paper:
% exogenous: usigmar usigmatb ur utb ug 
varlist_FV_exo_states={'u_sigma_r';'u_sigma_tb';'u_r';'u_tb';'u_x'}; 
for ii=1: length(varlist_FV_exo_states)
    FV_exo_order(ii)=strmatch(varlist_FV_exo_states(ii),M_.exo_names,'exact');
end

%followed by endogenous: sigmarlag sigmatblag erlag etblag glag investlag debt capital Lambda        
varlist_FV_endo_states={'sigma_r';'sigma_tb';'eps_r';'eps_tb';'X';'I';'D';'K'};
for ii=1: length(varlist_FV_endo_states)
    FV_endo_state_order(ii)=strmatch(varlist_FV_endo_states(ii),M_.endo_names(state_vars,:),'exact');
end

%% First order
gx_dyn(:,1:nu)=dr.ghu(FV_endo_order,FV_exo_order);
gx_dyn(:,nu+1:nu+nx)=dr.ghx(FV_endo_order,FV_endo_state_order);
if max(max(abs(gx(:,1:end-1)-gx_dyn))) > 1e-9
    max(max(abs(gx(:,1:end-1)-gx_dyn)))
    error('First order wrong')
else max(max(abs(gx(:,1:end-1)-gx_dyn)));
end

%% Second order
gxx_dyn=zeros(size(gxx));
for endo_iter_1=1:nx
    for endo_iter_2=1:nx
        gxx_dyn(nu+endo_iter_1,nu+endo_iter_2,:)=dr.ghxx(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nx+FV_endo_state_order(endo_iter_2));
    end
end

for exo_iter_1=1:nu
    for exo_iter_2=1:nu
        gxx_dyn(exo_iter_1,exo_iter_2,:)=dr.ghuu(FV_endo_order,(FV_exo_order(exo_iter_1)-1)*nu+FV_exo_order(exo_iter_2));
    end
end

for endo_iter_1=1:nx
    for exo_iter_1=1:nu
        gxx_dyn(nu+endo_iter_1,exo_iter_1,:)=dr.ghxu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nu+FV_exo_order(exo_iter_1));
        gxx_dyn(exo_iter_1,nu+endo_iter_1,:)=dr.ghxu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nu+FV_exo_order(exo_iter_1));
    end
end

% deal with Lambda, the perturbation parameter
for ii=1:length(FV_endo_order)
    gxx_dyn(14,14,:)=dr.ghs2(FV_endo_order);
end

if max(max(max(abs(gxx-gxx_dyn)))) > 1e-9
    max(max(max(abs(gxx-gxx_dyn))))
    error('Second order wrong')
else max(max(max(abs(gxx-gxx_dyn))));
end


%% third order
gxxx_dyn=zeros(size(gxxx));
for endo_iter_1=1:nx
    for endo_iter_2=1:nx
         for endo_iter_3=1:nx
            gxxx_dyn(nu+endo_iter_1,nu+endo_iter_2,nu+endo_iter_3,:)=dr.ghxxx(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nx*nx+(FV_endo_state_order(endo_iter_2)-1)*nx+FV_endo_state_order(endo_iter_3));
         end
    end
end

for exo_iter_1=1:nu
    for exo_iter_2=1:nu
        for exo_iter_3=1:nu
            gxxx_dyn(exo_iter_1,exo_iter_2,exo_iter_3,:)=dr.ghuuu(FV_endo_order,(FV_exo_order(exo_iter_1)-1)*nu*nu+(FV_exo_order(exo_iter_2)-1)*nu+FV_exo_order(exo_iter_3));
        end
    end
end

for endo_iter_1=1:nx
    for endo_iter_2=1:nx
         for exo_iter=1:nu
            gxxx_dyn(nu+endo_iter_1,nu+endo_iter_2,exo_iter,:)=dr.ghxxu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nx*nu+(FV_endo_state_order(endo_iter_2)-1)*nu+FV_exo_order(exo_iter));
            gxxx_dyn(exo_iter,nu+endo_iter_2,nu+endo_iter_1,:)=dr.ghxxu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nx*nu+(FV_endo_state_order(endo_iter_2)-1)*nu+FV_exo_order(exo_iter));
            gxxx_dyn(nu+endo_iter_1,exo_iter,nu+endo_iter_2,:)=dr.ghxxu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nx*nu+(FV_endo_state_order(endo_iter_2)-1)*nu+FV_exo_order(exo_iter));      
         end
    end
end

for endo_iter_1=1:nx
    for exo_iter_1=1:nu
         for exo_iter_2=1:nu
            gxxx_dyn(nu+endo_iter_1,exo_iter_1,exo_iter_2,:)=dr.ghxuu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nu*nu+(FV_exo_order(exo_iter_1)-1)*nu+FV_exo_order(exo_iter_2));
            gxxx_dyn(exo_iter_1,nu+endo_iter_1,exo_iter_2,:)=dr.ghxuu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nu*nu+(FV_exo_order(exo_iter_1)-1)*nu+FV_exo_order(exo_iter_2));
            gxxx_dyn(exo_iter_1,exo_iter_2,nu+endo_iter_1,:)=dr.ghxuu(FV_endo_order,(FV_endo_state_order(endo_iter_1)-1)*nu*nu+(FV_exo_order(exo_iter_1)-1)*nu+FV_exo_order(exo_iter_2));
         end
    end
end

% deal with Lambda, the perturbation parameter
gxxx_dyn(14,14,1:5,:)=oo_.dr.ghuss(FV_endo_order,FV_exo_order)';
gxxx_dyn(14,1:5,14,:)=oo_.dr.ghuss(FV_endo_order,FV_exo_order)';
gxxx_dyn(1:5,14,14,:)=oo_.dr.ghuss(FV_endo_order,FV_exo_order)';

gxxx_dyn(14,14,nu+1:13,:)=oo_.dr.ghxss(FV_endo_order,FV_endo_state_order)';
gxxx_dyn(14,14,14,:)=0;
gxxx_dyn(14,nu+1:13,14,:)=oo_.dr.ghxss(FV_endo_order,FV_endo_state_order)';
gxxx_dyn(nu+1:13,14,14,:)=oo_.dr.ghxss(FV_endo_order,FV_endo_state_order)';

if max(max(max(max(abs(gxxx-gxxx_dyn))))) > 1e-9
    max(max(max(max(abs(gxxx-gxxx_dyn)))))
    error('Third order wrong')
else max(max(max(max(abs(gxxx-gxxx_dyn)))))
end

