//  Based on Luca Marchiori/Olivier Pierrard(2012) LOLA 2.0: Luxembourg OverLapping generation model for policy Analysis
// Involves a call to solve_one_boundary.m that is tested here

load lola_data.mat

% ====================================================
% declarations var -- varexo -- para   
% ====================================================

@#define nbr_work_generations=9
@#define nbr_early_generations=2
@#define nbr_generations=16

parameters
length_period age_early;

length_period=5;
age_early=55;

@#define wt=[1]

@#define wg=0:nbr_work_generations-1
@#define ag=0:nbr_generations-1
@#define fwg=0:nbr_work_generations-nbr_early_generations-1
@#define nbwg=1:nbr_work_generations-1
@#define nbg=1:nbr_generations-1
@#define rg=nbr_work_generations:nbr_generations-1
@#define erg=nbr_work_generations-nbr_early_generations:nbr_work_generations-1
@#define endg=[nbr_generations-1]
@#define endw=[nbr_work_generations-1]

@#for i in wg
var 
n@{i} u@{i} Omega@{i} w@{i} dWHN@{i} dWFN@{i} 
n@{i}_f u@{i}_f Omega@{i}_f w@{i}_f dWFN@{i}_f
i@{i} lambda@{i} i@{i}_f lambda@{i}_f eta@{i};
parameters
Du@{i} Dn@{i} h@{i} h@{i}_f chi@{i} eta@{i}b; 
varexo
eps_eta@{i};
@#endfor

@#for i in ag
var 
c@{i} s@{i} P@{i} P@{i}_f;
varexo
beta@{i} beta@{i}_f PD@{i};
@#endfor

@#for i in erg
var 
WE@{i} De_@{i};
parameters
De_@{i}b; 
varexo
eps_De_@{i};
@#endfor

var   
wb wb_f 
Omega Omega_f Omega_hf 
V M qq p 
N N_f 
Q RR H K Y gdp nx FH pi
ct st wshare rr 
gamma mc phii D DH DF X bs bsY P00_f

rhou rhoe rhol tauw tauc tauf tauk g 
TFP gh rrb
theta tau1 om1 om2 om2s Ds phijs

DepRatio DepRatio_n DepRatio_d ZARA Ptot Ptot_f sleep du de dl inA inB in 
NBR NBRY NBR2 tauw2 tauf2 tauc2
PensCorr_L PensCorr_F;

parameters 
rho phi delta alpha beta ann 
fc nu aa 

rhoub rhoeb rholb tauwb taucb taufb taukb gb  
TFPb ghb rrbb
thetab tau1b om1b om2b om2sb Dsb phijsb

NBRYb bsY_iss;

varexo
P00 P00_foP00

eps_rhol eps_tauw eps_tauf eps_tauc eps_tauk
eps_rhoe eps_rhou eps_TFP eps_gh eps_theta eps_g
eps_Ds eps_phijs eps_PensCorr_L eps_PensCorr_F;


% ============================================================
% initialization
% ============================================================

@#for i in wg
set_param_value('Du@{i}',Du@{i});
set_param_value('Dn@{i}',Dn@{i});
set_param_value('h@{i}',h@{i});
set_param_value('h@{i}_f',h@{i}_f);
set_param_value('chi@{i}',chi@{i});
set_param_value('eta@{i}b',eta@{i}b);
@#endfor

@#for i in erg
set_param_value('De_@{i}b',De_@{i}b);
@#endfor

set_param_value('rho',rho);
set_param_value('phi',phi);
set_param_value('delta',delta);
set_param_value('alpha',alpha);
set_param_value('beta',beta);
set_param_value('ann',ann);
set_param_value('fc',fc);
set_param_value('nu',nu);
set_param_value('aa',aa);

set_param_value('rhoub',rhoub);
set_param_value('rhoeb',rhoeb);
set_param_value('rholb',rholb);
set_param_value('tauwb',tauwb);
set_param_value('taucb',taucb);
set_param_value('taufb',taufb);
set_param_value('taukb',taukb);
set_param_value('gb',gb);

set_param_value('TFPb',TFPb);
set_param_value('ghb',ghb);
set_param_value('rrbb',rrbb);

set_param_value('thetab',thetab);
set_param_value('tau1b',tau1b);
set_param_value('om1b',om1b);
set_param_value('om2b',om2b);
set_param_value('om2sb',om2sb);
set_param_value('Dsb',Dsb);
set_param_value('phijsb',phijsb);

set_param_value('bsY_iss',bsY_iss);

NBRYb=NBR_iss/(phii_iss*gdp_iss);


% =======================================================
model(block);              
% ======================================================

%  Labor Market Variables in the home country
%  ------------------------------------------ 

@#for i in fwg
0=lambda@{i};
@#endfor

@#for i in wg
1=n@{i}+u@{i}+i@{i};
@#endfor   
                            
i0=lambda0;
@#for i in nbwg
i@{i}=lambda@{i-1}(-1)+lambda@{i}*(1-lambda@{i-1}(-1));
@#endfor 

P0=beta0*P00+PD0;
@#for i in nbg
P@{i}=beta@{i}*P@{i-1}(-1)+PD@{i};
@#endfor

Omega0=P0;
@#for i in nbwg
Omega@{i}=(1-lambda@{i})*( 1-lambda@{i-1}(-1)-(1-chi@{i})*n@{i-1}(-1))*P@{i};
@#endfor 

n0=p;
@#for i in nbwg
n@{i}=(1-lambda@{i})*((1-p)*(1-chi@{i})*n@{i-1}(-1)+p*(1-lambda@{i-1}(-1)));
@#endfor 

N=
@#for i in wg
+n@{i}*P@{i}
@#endfor
;

%  Labor Market Variables in the foreign country
%  --------------------------------------------- 

@#for i in wg
1=n@{i}_f+u@{i}_f+i@{i}_f;
@#endfor    
                            
i0_f=lambda0_f;
@#for i in nbwg
i@{i}_f=lambda@{i-1}_f(-1)+lambda@{i}_f*(1-lambda@{i-1}_f(-1));
@#endfor

% -----------  reproduction cross-border   --------------------

P00_f=P00_foP00*P00;

P0_f=beta0_f*P00_f;
@#for i in nbg
P@{i}_f=beta@{i}_f*P@{i-1}_f(-1);
@#endfor 

Omega0_f=P0_f;
@#for i in nbwg
Omega@{i}_f=(1-lambda@{i}_f)*(1-lambda@{i-1}_f(-1)-(1-chi@{i})*n@{i-1}_f(-1))*P@{i}_f;
@#endfor 

n0_f=p;
@#for i in nbwg
n@{i}_f=(1-lambda@{i}_f)*((1-p)*(1-chi@{i})*n@{i-1}_f(-1)+p*(1-lambda@{i-1}_f(-1)));
@#endfor

N_f=
@#for i in wg
+n@{i}_f*P@{i}_f
@#endfor
;

%  Matching
% ----------

Omega=
@#for i in wg
+Omega@{i}
@#endfor
;

Omega_f=
@#for i in wg
+Omega@{i}_f
@#endfor
;

Omega_hf=Omega+Omega_f;
 
M=V*Omega_hf/(V^nu+Omega_hf^nu)^(1/nu);

qq=M/V;
p=M/Omega_hf;

%  Flow Budget Constraints (no bequests)
% --------------------------------------

rhou*w0*u0+ (1-tauw)*w0*n0 = (1+tauc)*c0+s0;
@#for i in nbwg
(1+rr*(1-tauk))/(beta@{i})^ann*s@{i-1}(-1)/(1+gh)+rhou*w@{i}*u@{i}+rhoe*w@{i}*i@{i}+(1-tauw)*w@{i}*n@{i}=(1+tauc)*c@{i}+s@{i};
@#endfor
@#for i in rg
(1+rr*(1-tauk))/(beta@{i})^ann*s@{i-1}(-1)/(1+gh)+rhol*wb=(1+tauc)*c@{i}+s@{i};
@#endfor
@#for i in endg
s@{i}=0;
@#endfor

wb=
@#for i in wg
+w@{i}/@{nbr_work_generations}
@#endfor
;

%  Euler Conditions
% ------------------

@#for i in nbg
1/(1+tauc)/c@{i-1}=beta*(1+rr(+1)*(1-tauk(+1)))/(1+tauc(+1))/c@{i}(+1)*(beta@{i})^(1-ann)/(1+gh);
@#endfor


%  Optimal Participation Rates (Early Retirement)
%  ----------------------------------------------

@#for i in erg
WE@{i} = 0;
@#endfor 

@#for i in erg
@#if i in endw
WE@{i} = ((De_@{i}b*i@{i}^(phi-1))+Du@{i}+(rhoe-rhou)*w@{i}/(1+tauc)/c@{i})*(1-i@{i})-((1-tauw-rhou)*w@{i}/(1+tauc)/c@{i}-(Dn@{i}-Du@{i}))*n@{i};
@#else  
WE@{i} = ((De_@{i}b*i@{i}^(phi-1))+Du@{i}+(rhoe-rhou)*w@{i}/(1+tauc)/c@{i})*(1-i@{i})-((1-tauw-rhou)*w@{i}/(1+tauc)/c@{i}-(Dn@{i}-Du@{i}))*n@{i}+ beta*beta@{i+1}*WE@{i+1};
@#endif
@#endfor  

%  Household Surplus
% -------------------

@#for i in wg
@#if i in endw
dWHN@{i} = w@{i}*(1-tauw-rhou)/(1+tauc)-(Dn@{i}-Du@{i})*c@{i};
@#else  
dWHN@{i} = w@{i}*(1-tauw-rhou)/(1+tauc)-(Dn@{i}-Du@{i})*c@{i} + beta*beta@{i+1}*c@{i}/c@{i+1}(+1)*dWHN@{i+1}(+1)*(1-p(+1))*(1-chi@{i+1})*(1-lambda@{i+1}(+1));
@#endif
@#endfor    

%  Foreign household
%  ------------------------
%  participation and wages
% ........................

@#for i in wg
lambda@{i}=lambda@{i}_f;
w@{i}_f=w@{i}; 
@#endfor

wb_f=
@#for i in wg
+w@{i}_f/@{nbr_work_generations}
@#endfor
;

%  Firm's Behavior
% -------------------

H=
@#for i in wg
+h@{i}*n@{i}*P@{i}+h@{i}_f*n@{i}_f*P@{i}_f
@#endfor
;

wshare=(1+tauf)*(
@#for i in wg
+w@{i}*n@{i}*P@{i}+w@{i}_f*n@{i}_f*P@{i}_f
@#endfor
)/(phii*gdp);

Y=TFP*H^(1-alpha)*(K)^alpha;
gdp= TFP*H^(1-alpha)*(K)^alpha-aa*V/phii-fc/phii;
pi = phii*gdp - wshare*phii*gdp - (rr+delta)*K(-1)/(1+gh);
(rr(+1)+delta)/(1+rr(+1)*(1-tauk(+1))) = mc*TFP*alpha*(H/(K))^(1-alpha);
FH=TFP*(1-alpha)*((K)/H)^alpha;

RR=1+rr*(1-tauk);
rr=rrb+tau1*(exp(bsY_iss-bsY)-1);

%  Firm's Surplus
%  ---------------------

@#for i in wg
@#if i in endw
dWFN@{i} = h@{i}*mc*FH-(1+tauf)*w@{i};
dWFN@{i}_f = h@{i}_f*mc*FH-(1+tauf)*w@{i}_f;
@#else  
dWFN@{i} = h@{i}*mc*FH-(1+tauf)*w@{i} + beta@{i+1}/RR(+1)*dWFN@{i+1}(+1)*(1-chi@{i+1})*(1-lambda@{i+1}(+1))*(1+gh);
dWFN@{i}_f = h@{i}_f*mc*FH-(1+tauf)*w@{i}_f + beta@{i+1}_f/RR(+1)*dWFN@{i+1}_f(+1)*(1-chi@{i+1})*(1-lambda@{i+1}_f(+1))*(1+gh);
@#endif
@#endfor

%  Free Entry Condition
%  ---------------------

aa=qq/Omega_hf*(
@#for i in wg
+Omega@{i}*dWFN@{i}+Omega@{i}_f*dWFN@{i}_f
@#endfor
);

%  Wage Determination (Rent Sharing)
%  -----------------------------------

@#for i in wg
(1-eta@{i})*dWHN@{i}  =  eta@{i}*((1-tauw)/(1+tauf)/(1+tauc))*dWFN@{i};
@#endfor

%  Equilibrium Conditions
%  ----------------------

ct=
@#for i in ag
+c@{i}*P@{i}
@#endfor
;

st=
@#for i in ag
+s@{i}*P@{i}
@#endfor
;

%  Non-Arbitrage condition (physical capital-shares)
% .........................

Q(+1)+pi(+1)=(1+rr(+1))*Q/(1+gh);

%  New Open Economy Macroeconomics (NOEM)
%  ---------------------------------------

phii=mc/theta;
D= ct + K-(1-delta)*K(-1)/(1+gh)  + g*gdp*phii + fc+aa*V;
DH=(1/om1*phii)^(1/(rho-1))*D;
X=(1/om2s*phii/gamma)^(1/(rho-1))*Ds;
DF=(1/om2*gamma*phijs)^(1/(rho-1))*D;
nx=phii*X-phijs*gamma*DF;     
bsY=bs/(phii*gdp);
Y=DH+X;
phii*gdp=ct + K-(1-delta)*K(-1)/(1+gh) + g*gdp*phii+nx;

st=K+Q +bs ;

%  Policies
%  ----------

rhou=rhoub*eps_rhou;
rhoe=rhoeb*eps_rhoe;
rhol=rholb*eps_rhol;
g=gb*eps_g;

@#for i in erg
De_@{i}=De_@{i}b*eps_De_@{i};
@#endfor

TFP=TFPb*eps_TFP;
gh=ghb*eps_gh;

@#for i in wg
eta@{i}=eta@{i}b*eps_eta@{i};
@#endfor

rrb=rrbb;

theta=thetab*eps_theta;
tau1=tau1b;
om1=om1b;
om2=om2b;
om2s=om2sb;
Ds=Dsb*eps_Ds;
phijs=phijsb*eps_phijs;


% ----------- RefDR scenario 

DepRatio_n=
@#for i in rg
+P@{i}
@#endfor
;
DepRatio_d=
@#for i in wg
+P@{i}
@#endfor
;

DepRatio=DepRatio_n/DepRatio_d;

ZARA=age_early+length_period*(
@#for i in erg
+1-i@{i}
@#endfor  
);

% ----------- WGEM 

Ptot=
@#for i in ag
+P@{i}
@#endfor
;
Ptot_f=
@#for i in ag
+P@{i}_f
@#endfor
;

sleep=(1+rr)*(
@#for i in nbg
+1/beta@{i}*(1-1/beta@{i}^(ann-1))*s@{i-1}(-1)*P@{i}
@#endfor  
)/(1+gh);

du=rhou*(
@#for i in wg
+w@{i}*u@{i}*P@{i}
@#endfor  
);

de=rhoe*(
@#for i in erg
+w@{i}*i@{i}*P@{i}+w@{i}_f*i@{i}_f*P@{i}_f
@#endfor  
);

dl=
@#for i in rg
+rhol*wb*PensCorr_L*P@{i}+rhol*wb_f*(N_f/(N+N_f))*PensCorr_F*P@{i}_f
@#endfor  
;

PensCorr_L=eps_PensCorr_L; 
PensCorr_F=eps_PensCorr_F;   

inA=(tauw+tauf)*(
@#for i in wg
+n@{i}*w@{i}*P@{i}+n@{i}_f*w@{i}_f*P@{i}_f
@#endfor  
);

inB=tauk*rr*(
@#for i in nbg
+1/beta@{i}^ann*s@{i-1}(-1)*P@{i}
@#endfor  
)/(1+gh);

in=tauc*ct+inA+inB+sleep;
NBR=g*phii*gdp+(du+de+dl)-(in);
NBRY=NBR/(phii*gdp);

% ----------- WGEM Adjustment variable ---------------

tauf2=tauf;
tauw2=tauw;
NBR2=NBR;
tauc2=tauc;
tauf=taufb*eps_tauf;
%----- WGEM: adjustment through tauc
tauc=taucb*eps_tauc;
%----- WGEM: adjustment through tauk
tauk=taukb*eps_tauk;
%----- WGEM: adjustment through tauw
tauw=tauwb*eps_tauw;

end;

%==================================
initval;
%==================================

@#for i in wg
n@{i}=n@{i}_iss;
n@{i}_f=n@{i}_f_iss;
u@{i}=u@{i}_iss;
u@{i}_f=u@{i}_f_iss;
Omega@{i}=Omega@{i}_iss;
Omega@{i}_f=Omega@{i}_f_iss;
w@{i}=w@{i}_iss;
w@{i}_f=w@{i}_f_iss;
dWHN@{i}=dWHN@{i}_iss;
dWFN@{i}=dWFN@{i}_iss;
dWFN@{i}_f=dWFN@{i}_f_iss;
eps_eta@{i}=eps_eta@{i}_iss;
eta@{i}=eta@{i}b*eps_eta@{i};
@#endfor

@#for i in erg
i@{i}=i@{i}_iss;
lambda@{i}=lambda@{i}_iss;
i@{i}_f=i@{i}_f_iss;
lambda@{i}_f=lambda@{i}_f_iss;
WE@{i}=0;
eps_De_@{i}=eps_De_@{i}_iss;
De_@{i}=De_@{i}b*eps_De_@{i}_iss;
@#endfor

@#for i in fwg
i@{i}=0;
lambda@{i}=0;
i@{i}_f=0;
lambda@{i}_f=0;
@#endfor

@#for i in ag
@#if i in endg
s@{i}=0;
@#else
s@{i}=s@{i}_iss;
@#endif
@#endfor

@#for i in ag
beta@{i}=beta@{i}_iss;
beta@{i}_f=beta@{i}_f_iss;
PD@{i}=PD@{i}_iss;
c@{i}=c@{i}_iss;
P@{i}=P@{i}_iss;
P@{i}_f=P@{i}_f_iss;
@#endfor

wb          =   wb_iss;
wb_f        =   wb_f_iss;
Omega       =   Omega_iss;
Omega_f     =   Omega_f_iss;
Omega_hf    =   Omega_hf_iss;
V        	=	V_iss	;
M        	=	M_iss	;
qq       	=	qq_iss	;
p        	=	p_iss	;
N        	=	N_iss	;
N_f      	=	N_f_iss	;
Q        	=	Q_iss	;
RR       	=	RR_iss	;
H        	=	H_iss	;
K        	=	K_iss	;
Y        	=	Y_iss	;
gdp      	=	gdp_iss	;
nx       	=	nx_iss	;
FH       	=	FH_iss	;
pi       	=	pi_iss	;
ct       	=	ct_iss	;
st       	=	st_iss	;
wshare   	=	wshare_iss	;
rr       	=	rr_iss	;

gamma    	=	gamma_iss	;
mc       	=	mc_iss	;
phii     	=	phii_iss	;
D        	=	D_iss	;
DH       	=	DH_iss	;
DF       	=	DF_iss	;
X        	=	X_iss	;
bs       	=	bs_iss	;
bsY      	=	bsY_iss	;
P00_f       =   P00_f_iss;

eps_rhol=eps_rhol_iss;
eps_rhoe=eps_rhoe_iss;
eps_rhou=eps_rhou_iss;
rhou=rhoub*eps_rhou_iss;
rhoe=rhoeb*eps_rhoe_iss;
rhol=rholb*eps_rhol_iss;
eps_tauc=eps_tauc_iss;
eps_tauk=eps_tauk_iss;
eps_tauw=eps_tauw_iss;
eps_tauf=eps_tauf_iss;
tauw=tauwb*eps_tauw;
tauc=taucb*eps_tauc;
tauf=taufb*eps_tauf;
tauk=taukb*eps_tauk;
eps_theta=eps_theta_iss;
eps_gh=eps_gh_iss;
eps_TFP=eps_TFP_iss;
eps_g=eps_g_iss;
g=gb*eps_g_iss;
TFP=TFPb*eps_TFP_iss;
gh=ghb*eps_gh_iss;
theta=thetab*eps_theta_iss;

rrb=rrbb;
tau1=tau1b;
om1=om1b;
om2=om2b;
om2s=om2sb;

eps_Ds=1;
eps_phijs=1;
Ds=Dsb*eps_Ds;
phijs=phijsb*eps_phijs;

eps_PensCorr_F=eps_PensCorr_F_iss;
eps_PensCorr_L=eps_PensCorr_L_iss;
PensCorr_F=eps_PensCorr_F_iss;
PensCorr_L=eps_PensCorr_L_iss;

P00	        =	P00_iss	;
P00_foP00	=	P00_foP00_iss	;

DepRatio_n=
@#for i in rg
+P@{i}
@#endfor
;
DepRatio_d=
@#for i in wg
+P@{i}
@#endfor
;

DepRatio=DepRatio_n/DepRatio_d;

ZARA=age_early+length_period*(
@#for i in erg
+1-i@{i}
@#endfor  
);

Ptot=Ptot_iss;
Ptot_f=Ptot_f_iss;
sleep=sleep_iss;
du=du_iss;
de=de_iss;
dl=dl_iss;

inA=inA_iss;%(tauw+tauf)*(
%@#for i in wg
%+n@{i}*w@{i}*P@{i}+n@{i}_f*w@{i}_f*P@{i}_f
%@#endfor  
%);

inB=inB_iss;%tauk*rr*(
%@#for i in nbg
%+1/beta@{i}^ann*s@{i-1}*P@{i}
%@#endfor  
%)/(1+gh);

in=in_iss;%tauc*ct+inA+inB+sleep;
NBR=NBR_iss;%g*phii*gdp+(du+de+dl)-(in);
NBRY=NBR/(phii*gdp);
NBR2=NBR;
tauf2=tauf;
tauw2=tauw;
tauc2=tauc;

end;

%========================================================
% compute initial steady state and check eigenvalues 
%========================================================

resid;
steady(solve_algo=3);
check;

%========================================================
endval;
%========================================================

@#for i in wg
n@{i}=n@{i}_fss;
n@{i}_f=n@{i}_f_fss;
u@{i}=u@{i}_fss;
u@{i}_f=u@{i}_f_fss;
Omega@{i}=Omega@{i}_fss;
Omega@{i}_f=Omega@{i}_f_fss;
w@{i}=w@{i}_fss;
w@{i}_f=w@{i}_f_fss;
dWHN@{i}=dWHN@{i}_fss;
dWFN@{i}=dWFN@{i}_fss;
dWFN@{i}_f=dWFN@{i}_f_fss;
eps_eta@{i}=eps_eta@{i}_fss;
eta@{i}=eta@{i}b*eps_eta@{i};
@#endfor

@#for i in erg
i@{i}=i@{i}_fss;
lambda@{i}=lambda@{i}_fss;
i@{i}_f=i@{i}_f_fss;
lambda@{i}_f=lambda@{i}_f_fss;
WE@{i}=0;
eps_De_@{i}=eps_De_@{i}_fss;
De_@{i}=De_@{i}b*eps_De_@{i}_fss;
@#endfor

@#for i in fwg
i@{i}=0;
lambda@{i}=0;
i@{i}_f=0;
lambda@{i}_f=0;
@#endfor

@#for i in ag
@#if i in endg
s@{i}=0;
@#else
s@{i}=s@{i}_fss;
@#endif
@#endfor

@#for i in ag
beta@{i}=beta@{i}_fss;
beta@{i}_f=beta@{i}_f_fss;
PD@{i}=PD@{i}_fss;
c@{i}=c@{i}_fss;
P@{i}=P@{i}_fss;
P@{i}_f=P@{i}_f_fss;
@#endfor

wb          =   wb_fss;
wb_f        =   wb_f_fss;
Omega       =   Omega_fss;
Omega_f     =   Omega_f_fss;
Omega_hf    =   Omega_hf_fss;
V        	=	V_fss	;
M        	=	M_fss	;
qq       	=	qq_fss	;
p        	=	p_fss	;
N        	=	N_fss	;
N_f      	=	N_f_fss	;
Q        	=	Q_fss	;
RR       	=	RR_fss	;
H        	=	H_fss	;
K        	=	K_fss	;
Y        	=	Y_fss	;
gdp      	=	gdp_fss	;
nx       	=	nx_fss	;
FH       	=	FH_fss	;
pi       	=	pi_fss	;
ct       	=	ct_fss	;
st       	=	st_fss	;
wshare   	=	wshare_fss	;
rr       	=	rr_fss	;

gamma    	=	gamma_fss	;
mc       	=	mc_fss	;
phii     	=	phii_fss	;
D        	=	D_fss	;
DH       	=	DH_fss	;
DF       	=	DF_fss	;
X        	=	X_fss	;
bs       	=	bs_fss	;
bsY      	=	bsY_fss	;
P00_f       =   P00_f_fss;

eps_rhol=eps_rhol_fss;
eps_rhoe=eps_rhoe_fss;
eps_rhou=eps_rhou_fss;
rhou=rhoub*eps_rhou_fss;
rhoe=rhoeb*eps_rhoe_fss;
rhol=rholb*eps_rhol_fss;
eps_tauc=eps_tauc_fss;
eps_tauk=eps_tauk_fss;
eps_tauw=eps_tauw_fss;
eps_tauf=eps_tauf_fss;
tauw=tauwb*eps_tauw;
tauc=taucb*eps_tauc;
tauf=taufb*eps_tauf;
tauk=taukb*eps_tauk;
eps_theta=eps_theta_fss;
eps_gh=eps_gh_fss;
eps_TFP=eps_TFP_fss;
eps_g=eps_g_fss;
g=gb*eps_g_fss;
TFP=TFPb*eps_TFP_fss;
gh=ghb*eps_gh_fss;
theta=thetab*eps_theta_fss;

rrb=rrbb;
tau1=tau1b;
om1=om1b;
om2=om2b;
om2s=om2sb;

eps_Ds=1;
eps_phijs=1;
Ds=Dsb*eps_Ds;
phijs=phijsb*eps_phijs;

eps_PensCorr_F=eps_PensCorr_F_fss;
eps_PensCorr_L=eps_PensCorr_L_fss;
PensCorr_F=eps_PensCorr_F_fss;
PensCorr_L=eps_PensCorr_L_fss;

P00	        =	P00_fss	;
P00_foP00	=	P00_foP00_fss	;

DepRatio_n=
@#for i in rg
+P@{i}
@#endfor
;
DepRatio_d=
@#for i in wg
+P@{i}
@#endfor
;

DepRatio=DepRatio_n/DepRatio_d;

ZARA=age_early+length_period*(
@#for i in erg
+1-i@{i}
@#endfor  
);

Ptot=
@#for i in ag
+P@{i}
@#endfor
;
Ptot_f=
@#for i in ag
+P@{i}_f
@#endfor
;

sleep=(1+rr)*(
@#for i in nbg
+1/beta@{i}*(1-1/beta@{i}^(ann-1))*s@{i-1}*P@{i}
@#endfor  
)/(1+gh);

du=rhou*(
@#for i in wg
+w@{i}*u@{i}*P@{i}
@#endfor  
);

de=rhoe*(
@#for i in erg
+w@{i}*i@{i}*P@{i}+w@{i}_f*i@{i}_f*P@{i}_f
@#endfor  
);

dl=
@#for i in rg
+rhol*wb*PensCorr_L*P@{i}+rhol*wb_f*(N_f/(N+N_f))*PensCorr_F*P@{i}_f
@#endfor  
;

inA=(tauw+tauf)*(
@#for i in wg
+n@{i}*w@{i}*P@{i}+n@{i}_f*w@{i}_f*P@{i}_f
@#endfor  
);

inB=tauk*rr*(
@#for i in nbg
+1/beta@{i}^ann*s@{i-1}*P@{i}
@#endfor  
)/(1+gh);

in=tauc*ct+inA+inB+sleep;
NBR=g*phii*gdp+(du+de+dl)-(in);
NBRY=NBR/(phii*gdp);
NBR2=NBR;
tauf2=tauf;
tauw2=tauw;
tauc2=tauc;

end;

%========================================================
% compute final steady state and check eigenvalues 
%========================================================

resid;
steady(solve_algo=3);
check;


% ===================================================
shocks;     
% ===================================================

var P00;																																																			
   periods 1:99;																																																			
   values (se_P00);																																																			
 
@#for i in nbg
var beta@{i};																																																			
   periods  1:99;																																																			
   values (se_beta@{i});
var beta@{i}_f;																																																			
   periods  1:99;																																																			
   values (se_beta@{i}_f);
var PD@{i};																																																			
   periods  1:99;																																																			
   values (se_PD@{i});
@#endfor  
																																																			
var P00_foP00;																																																			
   periods 1:99;																																																			
   values (se_P00_foP00);	 																																																

var eps_g;																																																			
   periods 1:99;																																																			
   values (se_eps_g);	

var eps_PensCorr_F;																																																			
   periods 1:99;																																																			
   values (se_eps_PensCorr_F);	

var eps_PensCorr_L;																																																			
   periods 1:99;																																																			
   values (se_eps_PensCorr_L);	

end;

% *******************************************
% Numerical Simulation, Control Parameters
% *******************************************

simul(periods=125,maxit=100);



