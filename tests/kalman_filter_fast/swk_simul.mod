// country_nbr must be passed on the command line: dynare swk_0 -Dcountry_nbr=3
@#define countries = [1:country_nbr]

var 
@#for c in countries
MC_@{c} E_@{c} EF_@{c} R_KF_@{c} QF_@{c} CF_@{c} IF_@{c} YF_@{c} LF_@{c} PIEF_@{c} WF_@{c} RF_@{c} R_K_@{c} Q_@{c} C_@{c} I_@{c} Y_@{c} L_@{c} PIE_@{c} W_@{c} R_@{c} EE_A_@{c} PIE_BAR_@{c} EE_B_@{c} EE_G_@{c} EE_L_@{c} EE_I_@{c} KF_@{c} K_@{c}
@#endfor    
;
varexo 
@#for c in countries
E_A_@{c} E_B_@{c} E_G_@{c} E_L_@{c} E_I_@{c} ETA_R_@{c} E_PIE_BAR_@{c} ETA_Q_@{c} ETA_P_@{c} ETA_W_@{c}
@#endfor
;
 
parameters xi_e lambda_w alpha czcap beta phi_i tau sig_c hab ccs cinvs phi_y gamma_w xi_w gamma_p xi_p sig_l r_dpi r_pie r_dy r_y rho rho_a rho_pb rho_b rho_g rho_l rho_i  ;
alpha=.30;
beta=.99;
tau=0.025;
ccs=0.6;
cinvs=.22; 
lambda_w = 0.5;
phi_i= 6.771;
sig_c=   1.353; 
hab=    0.573;    
xi_w=   0.737;
sig_l=    2.400;
xi_p=   0.908;
xi_e= 0.599;
gamma_w=    0.763;
gamma_p=    0.469;
czcap=    0.169;
phi_y=    1.408;
r_pie=     1.684;
r_dpi=    0.14;
rho=      0.961;
r_y=      0.099;
r_dy=     0.159;
rho_a=    0.823;
rho_b=    0.855;
rho_g=    0.949;
rho_l=   0.889;
rho_i=   0.927;
rho_pb=  0.924;


model(linear); 
@#for c in countries
          CF_@{c} = (1/(1+hab))*(CF_@{c}(1)+hab*CF_@{c}(-1))-((1-hab)/((1+hab)*sig_c))*(RF_@{c}-PIEF_@{c}(1)-EE_B_@{c}) ;
	      0 =  alpha*R_KF_@{c}+(1-alpha)*WF_@{c} -EE_A_@{c} ;
          PIEF_@{c} = 0;
	      IF_@{c} = (1/(1+beta))* ((  IF_@{c}(-1) + beta*(IF_@{c}(1)))+(1/phi_i)*QF_@{c})+0*ETA_Q_@{c}+EE_I_@{c} ;
	      QF_@{c} = -(RF_@{c}-PIEF_@{c}(1))+(1-beta*(1-tau))*((1+czcap)/czcap)*R_KF_@{c}(1)+beta*(1-tau)*QF_@{c}(1) +0*EE_I_@{c} ;
          KF_@{c} =  (1-tau)*KF_@{c}(-1)+tau*IF_@{c}(-1) ;
	      YF_@{c} = (ccs*CF_@{c}+cinvs*IF_@{c})+EE_G_@{c}  ;
	      YF_@{c} = 1*phi_y*( alpha*KF_@{c}+alpha*(1/czcap)*R_KF_@{c}+(1-alpha)*LF_@{c}+EE_A_@{c} ) ;
	      WF_@{c} = (sig_c/(1-hab))*(CF_@{c}-hab*CF_@{c}(-1)) + sig_l*LF_@{c} - EE_L_@{c} ;
	      LF_@{c} = R_KF_@{c}*((1+czcap)/czcap)-WF_@{c}+KF_@{c} ;
          EF_@{c} = EF_@{c}(-1)+EF_@{c}(1)-EF_@{c}+(LF_@{c}-EF_@{c})*((1-xi_e)*(1-xi_e*beta)/(xi_e));
	      C_@{c} = (hab/(1+hab))*C_@{c}(-1)+(1/(1+hab))*C_@{c}(1)-((1-hab)/((1+hab)*sig_c))*(R_@{c}-PIE_@{c}(1)-EE_B_@{c}) ;
	      I_@{c} = (1/(1+beta))* ((  I_@{c}(-1) + beta*(I_@{c}(1)))+(1/phi_i)*Q_@{c} )+1*ETA_Q_@{c}+1*EE_I_@{c} ;
	      Q_@{c} = -(R_@{c}-PIE_@{c}(1))+(1-beta*(1-tau))*((1+czcap)/czcap)*R_K_@{c}(1)+beta*(1-tau)*Q_@{c}(1) +EE_I_@{c}*0+0*ETA_Q_@{c} ;
	      K_@{c} =  (1-tau)*K_@{c}(-1)+tau*I_@{c}(-1) ;
	      Y_@{c} = (ccs*C_@{c}+cinvs*I_@{c})+ EE_G_@{c}   ;
	      Y_@{c} = phi_y*( alpha*K_@{c}+alpha*(1/czcap)*R_K_@{c}+(1-alpha)*L_@{c} ) +phi_y*EE_A_@{c}  ;
	      PIE_@{c} = (1/(1+beta*gamma_p))*((beta)*(PIE_@{c}(1)) +(gamma_p)*(PIE_@{c}(-1))+((1-xi_p)*(1-beta*xi_p)/(xi_p))*(MC_@{c}))  + ETA_P_@{c} ; 
	      MC_@{c} = alpha*R_K_@{c}+(1-alpha)*W_@{c} -EE_A_@{c};
	      W_@{c} =  (1/(1+beta))*(beta*W_@{c}(+1)+W_@{c}(-1))+(beta/(1+beta))*(PIE_@{c}(+1))-((1+beta*gamma_w)/(1+beta))*(PIE_@{c})+(gamma_w/(1+beta))*(PIE_@{c}(-1))-(1/(1+beta))*(((1-beta*xi_w)*(1-xi_w))/(((1+(((1+lambda_w)*sig_l)/(lambda_w))))*xi_w))*(W_@{c}-sig_l*L_@{c}-(sig_c/(1-hab))*(C_@{c}-hab*C_@{c}(-1))+EE_L_@{c})+ETA_W_@{c};
	      L_@{c} = R_K_@{c}*((1+czcap)/czcap)-W_@{c}+K_@{c} ;
	      R_@{c} = r_dpi*(PIE_@{c}-PIE_@{c}(-1))+(1-rho)*(r_pie*(PIE_@{c}(-1)-PIE_BAR_@{c})+r_y*(Y_@{c}-YF_@{c}))+r_dy*(Y_@{c}-YF_@{c}-(Y_@{c}(-1)-YF_@{c}(-1)))+rho*(R_@{c}(-1)-PIE_BAR_@{c})+PIE_BAR_@{c}+ETA_R_@{c};
          E_@{c} = E_@{c}(-1)+E_@{c}(1)-E_@{c}+(L_@{c}-E_@{c})*((1-xi_e)*(1-xi_e*beta)/(xi_e));
          EE_A_@{c} = (rho_a)*EE_A_@{c}(-1)  + E_A_@{c};
	      PIE_BAR_@{c} = rho_pb*PIE_BAR_@{c}(-1)+ E_PIE_BAR_@{c} ;
	      EE_B_@{c} = rho_b*EE_B_@{c}(-1) + E_B_@{c} ;
	      EE_G_@{c} = rho_g*EE_G_@{c}(-1) + E_G_@{c} ;
	      EE_L_@{c} = rho_l*EE_L_@{c}(-1) + E_L_@{c} ;
	      EE_I_@{c} = rho_i*EE_I_@{c}(-1) + E_I_@{c} ;
@#endfor
end; 

 
shocks;
@#for c in countries
var E_A_@{c}; stderr 0.543;
var E_PIE_BAR_@{c}; stderr 0.072;
var E_B_@{c}; stderr 0.2694;
var E_G_@{c}; stderr 0.3052;
var E_L_@{c}; stderr 1.4575;
var E_I_@{c}; stderr 0.1318;
var ETA_R_@{c}; stderr 0.1363;
var ETA_Q_@{c}; stderr 0.4842;
var ETA_P_@{c}; stderr 0.1731;
var ETA_W_@{c}; stderr 0.2462;
@#endfor
end;

stoch_simul(order=1,periods=200,irf=0,noprint);

@#for c in countries
save pseudo_data.mat Y_@{c} C_@{c} I_@{c} E_@{c} PIE_@{c} W_@{c} R_@{c} -append;
@#endfor


