var y y_s R pie dq pie_s de A y_obs pie_obs R_obs;
varexo e_R e_q e_ys e_pies e_A;

parameters psi1 psi2 psi3 rho_R tau alpha rr k rho_q rho_A rho_ys rho_pies;

psi1 = 1.54;
psi2 = 0.25;
psi3 = 0.25;
rho_R = 0.5;
alpha = 0.3;
rr = 2.51;
k = 0.5;
tau = 0.5;
rho_q = 0.4;
rho_A = 0.2;
rho_ys = 0.9;
rho_pies = 0.7;


model(linear);
y = y(+1) - (tau +alpha*(2-alpha)*(1-tau))*(R-pie(+1))-alpha*(tau +alpha*(2-alpha)*(1-tau))*dq(+1) + alpha*(2-alpha)*((1-tau)/tau)*(y_s-y_s(+1))-A(+1);
pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
pie = de+(1-alpha)*dq+pie_s;
R = rho_R*R(-1)+(1-rho_R)*(psi1*pie+psi2*(y+alpha*(2-alpha)*((1-tau)/tau)*y_s)+psi3*de)+e_R;
dq = rho_q*dq(-1)+e_q;
y_s = rho_ys*y_s(-1)+e_ys;
pie_s = rho_pies*pie_s(-1)+e_pies;
A = rho_A*A(-1)+e_A;
y_obs = y-y(-1)+A;
pie_obs = 4*pie;
R_obs = 4*R;
end;

shocks;
var e_R = 1.25^2;
var e_q = 2.5^2;
var e_A = 1.89;
var e_ys = 1.89;
var e_pies = 1.89;
end;

varobs y_obs R_obs pie_obs dq de;

estimated_params;
psi1 , gamma_pdf,1.5,0.5;
psi2 , gamma_pdf,0.25,0.125;
psi3 , gamma_pdf,0.25,0.125;
rho_R ,beta_pdf,0.5,0.2;
alpha ,beta_pdf,0.3,0.1;
rr ,gamma_pdf,2.5,1;
k , gamma_pdf,0.5,0.25;
tau ,gamma_pdf,0.5,0.2;
rho_q ,beta_pdf,0.4,0.2;
rho_A ,beta_pdf,0.5,0.2;
rho_ys ,beta_pdf,0.8,0.1;
rho_pies,beta_pdf,0.7,0.15;
stderr e_R,inv_gamma_pdf,1.2533,0.6551;
stderr e_q,inv_gamma_pdf,2.5066,1.3103;
stderr e_A,inv_gamma_pdf,1.2533,0.6551;
stderr e_ys,inv_gamma_pdf,1.2533,0.6551;
stderr e_pies,inv_gamma_pdf,1.88,0.9827;
end;

options_.TeX=1;
estimation(datafile='../ls2003/data_ca1',first_obs=8,nobs=79,mh_nblocks=10,prefilter=1,mh_jscale=0.5,mh_replic=0);
close all

shock_groups(name=trade);
supply = e_A ;
trade = e_q ;
monetary = e_R ;
end;

shock_groups(name=row);
supply = e_A ;
'RoW shocks' = e_q e_ys e_pies  ;
monetary = e_R ;
end;
options_.initial_date=dates('1989Q4'); % date arbitrarily set for testing purposes
shock_decomposition(use_shock_groups=trade) y_obs R_obs pie_obs dq de;

// various tests for plot_shock_decompositions
// standard plot [using trade group defined before]
plot_shock_decomposition(M_,oo_,options_,var_list_);

// test datailed, custom name and yoy plots
options_.shock_decomp.detail_plot = 1;
options_.shock_decomp.fig_names = 'MR';
options_.shock_decomp.type='yoy';
plot_shock_decomposition(M_,oo_,options_,var_list_);

close all,


// testing realtime decomposition
// first compute realtime decompositions [pre-processor not yet available]
options_.shock_decomp.forecast=8;
options_.shock_decomp.save_realtime = [5:4:options_.nobs]; % store values useful for annualized vars
oo_ = realtime_shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);

options_.shock_decomp.detail_plot = 0;
options_.shock_decomp.type='';
options_.use_shock_groups='';
options_.shock_decomp.fig_names = '';


//realtime pooled
options_.shock_decomp.realtime=1; 
options_.shock_decomp.vintage=0; 
plot_shock_decomposition(M_,oo_,options_,var_list_);

//conditional pooled
options_.shock_decomp.realtime=2; 
plot_shock_decomposition(M_,oo_,options_,var_list_);
options_.shock_decomp.realtime=2; 
// conditional 8-step ahead decomposition, given 1989q4
options_.shock_decomp.detail_plot = 1;
options_.shock_decomp.vintage=29; 
plot_shock_decomposition(M_,oo_,options_,var_list_);

close all,

//forecast pooled
options_.shock_decomp.detail_plot = 0;
options_.shock_decomp.realtime=3; 
options_.shock_decomp.vintage=0; 
plot_shock_decomposition(M_,oo_,options_,var_list_);

// forecast 8-step ahead decomposition, given 1989q4
options_.shock_decomp.detail_plot = 1;
options_.shock_decomp.realtime=3; 
options_.shock_decomp.vintage=29; 
plot_shock_decomposition(M_,oo_,options_,var_list_);

close all,

// now I test annualized variables
options_.shock_decomp.realtime=0; 
options_.shock_decomp.detail_plot = 1;
options_.shock_decomp.interactive=0;
options_.use_shock_groups='';
options_.shock_decomp.type='aoa'; 
options_.shock_decomp.q2a=1; 
options_.shock_decomp.islog=1; 
plot_shock_decomposition(M_,oo_,options_,'y');

options_.shock_decomp.realtime=1; 
options_.shock_decomp.vintage=0; 
plot_shock_decomposition(M_,oo_,options_,'y');

options_.shock_decomp.vintage=29; 
options_.shock_decomp.realtime=1; 
plot_shock_decomposition(M_,oo_,options_,'y');
options_.shock_decomp.realtime=2; 
plot_shock_decomposition(M_,oo_,options_,'y');
options_.shock_decomp.realtime=3; 
plot_shock_decomposition(M_,oo_,options_,'y');

close all

//test uimenu for groups
options_.shock_decomp.realtime=0; 
options_.shock_decomp.detail_plot = 1;
options_.shock_decomp.interactive=1;
options_.use_shock_groups='row';
options_.shock_decomp.type='qoq'; 
plot_shock_decomposition(M_,oo_,options_,var_list_);

options_.shock_decomp.detail_plot = 1;
options_.shock_decomp.realtime=3; 
options_.shock_decomp.vintage=29; 
plot_shock_decomposition(M_,oo_,options_,var_list_);

collect_latex_files;
if system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])
    error('TeX-File did not compile.')
end 
