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
pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+k*alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
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

  
disp(' ');
disp('NOW I DO STABILITY MAPPING and prepare sample for Reduced form Mapping');
disp(' ');
disp('Press ENTER to continue'); pause(5);

dynare_sensitivity(redform=1,nodisplay, graph_format=fig); //create sample of reduced form coefficients
// NOTE: since namendo is empty by default, 
// this call does not perform the mapping of reduced form coefficient: just prepares the sample

disp(' ');
disp('ANALYSIS OF REDUCED FORM COEFFICIENTS');
disp(' ');
disp('Press ENTER to continue'); pause(5);

dynare_sensitivity(nodisplay, graph_format=fig,load_stab=1,  // loead previously generated sample analysed for stability
redform=1,  // do the reduced form mapping
threshold_redform=[-1 0],  // filter reduced form coefficients (default=[])
namendo=(pie,R),  // evaluate relationships for pie and R (namendo=(:) for all variables)
namexo=(e_R),     // evaluate relationships with exogenous e_R (use namexo=(:) for all shocks)
namlagendo=(R),   // evaluate relationships with lagged R (use namlagendo=(:) for all lagged endogenous)
stab=0 // don't repeat again the stability mapping
);


disp(' ');
disp('THE PREVIOUS TWO CALLS COULD BE DONE TOGETHER');
disp('BY USING THE COMBINED CALL');
disp(' ');
disp('dynare_sensitivity(redform=1,')
disp('threshold_redform=[-1 0], namendo=(pie,R), namexo=(e_R), namlagendo=(R));')
disp(' ');
disp('Press ENTER to continue'); pause(5);
//dynare_sensitivity(
//redform=1, //create sample of reduced form coefficients
//logtrans_redform=1,  // estimate log-transformed reduced form coefficients (default=0)
//namendo=(pie,R),  // evaluate relationships for pie and R (namendo=(:) for all variables)
//namexo=(e_R),     // evaluate relationships with exogenous e_R (use namexo=(:) for all shocks)
//namlagendo=(R)   // evaluate relationships with lagged R (use namlagendo=(:) for all lagged endogenous)
//);



disp(' ');
disp('MC FILTERING(rmse=1), TO MAP THE FIT FROM PRIORS');
disp('Press ENTER to continue'); pause(5);

dynare_sensitivity(nodisplay, graph_format=fig,datafile=data_ca1,first_obs=8,nobs=79,prefilter=1, // also presample=2,loglinear, are admissible
load_stab=1,     // load prior sample
istart_rmse=2,   //start computing rmse from second observation (i.e. rmse does not inlude initial big error)
stab=0,          // don't  plot again stability analysis results
rmse=1          // do rmse analysis
);

disp(' ');
disp('THE PREVIOUS THREE CALLS COULD BE DONE TOGETHER');
disp('BY USING THE COMBINED CALL');
disp(' ');
disp('dynare_sensitivity(redform=1,')
disp('logtrans_redform=1, namendo=(pie,R), namexo=(e_R), namlagendo=(R),')   
disp('datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,')
disp('istart_rmse=2, rmse=1);')
disp(' ');
disp('Press ENTER to continue'); pause(5);
//dynare_sensitivity(
//redform=1, //create sample of reduced form coefficients
//logtrans_redform=1,  // estimate log-transformed reduced form coefficients (default=0)
//namendo=(pie,R),  // evaluate relationships for pie and R (namendo=(:) for all variables)
//namexo=(e_R),     // evaluate relationships with exogenous e_R (use namexo=(:) for all shocks)
//namlagendo=(R),   // evaluate relationships with lagged R (use namlagendo=(:) for all lagged endogenous)
//datafile=data_ca1,first_obs=8,nobs=79,prefilter=1, 
//istart_rmse=2,   //start computing rmse from second observation (i.e. rmse does not inlude initial big error)
//rmse=1,          // do rmse analysis
//);



disp(' ');
disp('I ESTIMATE THE MODEL');
disp(' ');
disp('Press ENTER to continue'); pause(5);

// run this to generate posterior mode and Metropolis files if not yet done
estimation(datafile=data_ca1,first_obs=8,nobs=79,mh_nblocks=2,
   prefilter=1,mh_jscale=0.5,mh_replic=5000, mode_compute=4, mh_drop=0.6, nodisplay, graph_format=fig,
   bayesian_irf, filtered_vars, smoother) y_obs R_obs pie_obs dq de;


// run this to produce posterior samples of filtered, smoothed and irf variables, if not yet done
//estimation(datafile=data_ca1,first_obs=8,nobs=79,mh_nblocks=2,prefilter=1,mh_jscale=0.3,
//          mh_replic=0, mode_file=ls2003_mode, mode_compute=0, load_mh_file, bayesian_irf,
//		  filtered_vars, smoother, mh_drop=0.6);

disp(' ');
disp('WE DO STABILITY MAPPING AGAIN, EXPLORING AROUND THE MODE BY A WIDTH +- 0.2 times the mode');
disp('THEN WE DO RMSE analysis in this neighborhood');
disp('Typical for ML estimation, also feasible for posterior mode');
disp(' ');
disp('Press ENTER to continue'); pause(5);

dynare_sensitivity(nodisplay, graph_format=fig,pprior=0,Nsam=2048,neighborhood_width=0.2,
mode_file=ls2003_mode,  // specifies the mode file where the mode and Hessian are stored
datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
rmse=1);

disp(' ');
disp('WE DO STABILITY MAPPING AGAIN, BUT FOR MULTIVARIATE SAMPLE AT THE POSTERIOR MODE (or ML) and Hessian (pprior=0 & ppost=0)');
disp('Typical for ML estimation, also feasible for posterior mode');
disp(' ');
disp('Press ENTER to continue'); pause(5);

dynare_sensitivity(nodisplay, graph_format=fig,pprior=0,Nsam=2048,
mode_file=ls2003_mode  // specifies the mode file where the mode and Hessian are stored
);


disp(' ');
disp('RMSE ANALYSIS FOR MULTIVARIATE SAMPLE AT THE POSTERIOR MODE');
disp(' ');
disp('Press ENTER to continue'); pause(5);
dynare_sensitivity(nodisplay, graph_format=fig,mode_file=ls2003_mode,
datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
pprior=0,
stab=0,
rmse=1,
alpha2_rmse=0, // no correlation analysis
alpha_rmse=0  // no Smirnov sensitivity analysis
);

disp(' ');
disp('THE LAST TWO CALLS COULD BE DONE TOGETHER');
disp('BY USING THE COMBINED CALL');
disp(' ');
disp('dynare_sensitivity(pprior=0,Nsam=2048,alpha2_stab=0.4,mode_file=ls2003_mode,')
disp('datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,')
disp('rmse=1, alpha2_rmse=0, alpha_rmse=0);')
disp(' ');
disp('Press ENTER to continue'); pause(5);
//dynare_sensitivity(pprior=0,Nsam=2048,alpha2_stab=0.4,mode_file=ls2003_mode,
//datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
//rmse=1,
//alpha2_rmse=0, // no correlation analysis
//alpha_rmse=0  // no Smirnov sensitivity analysis
//);

disp(' ');
disp('RMSE ANALYSIS FOR POSTERIOR MCMC sample (ppost=1)');
disp('Needs a call to dynare_estimation to load all MH environment');
disp('Press ENTER to continue'); pause(5);
//estimation(datafile=data_ca1,first_obs=8,nobs=79,mh_nblocks=2, mode_file=ls2003_mode, load_mh_file,
//  prefilter=1,mh_jscale=0.5,mh_replic=0, mode_compute=0, mh_drop=0.6);

dynare_sensitivity(nodisplay, graph_format=fig,stab=0, // no need for stability analysis since the posterior sample is surely OK
datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
rmse=1,ppost=1);

