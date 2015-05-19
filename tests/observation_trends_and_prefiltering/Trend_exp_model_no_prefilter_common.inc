var Y_obs P_obs junk1 junk2;
varexo e_y e_p eps_junk;

parameters rho_y rho_p  g_y g_p const_y const_p sigma_y sigma_p;

rho_y=0.5;
rho_p=0.5;
g_y=0.0001;
g_p=-0.0001;
const_y=2;
const_p=2;
sigma_y=0.001;
sigma_p=0.001;

model;
Y_obs = exp(const_y)^(1-rho_y)*Y_obs(-1)^rho_y*exp(sigma_y*e_y);
P_obs = exp(const_p)^(1-rho_p)*P_obs(-1)^rho_p*exp(sigma_p*e_p);
junk1 = (junk1(+1))^0.9;
junk2 = (junk2(-1))^0.9*exp(eps_junk);
end;

steady_state_model;
Y_obs = exp(const_y);
P_obs = exp(const_p);
junk1=1;
junk2=1;
end;

shocks;
var e_p; stderr 1;
var e_y; stderr 1;
var eps_junk; stderr 1;
end;

steady(nocheck);
check;

estimated_params;
const_y, normal_pdf, 2, 1;
const_p, normal_pdf, 2, 1;
g_y, normal_pdf, 0.0001, 1;
g_p, normal_pdf, -0.0001, 1;
rho_y, normal_pdf, 0.5, 1;
rho_p, normal_pdf, 0.5, 1;
sigma_y, inv_gamma_pdf, 0.001, inf;
sigma_p, inv_gamma_pdf, 0.001, inf;
end;

estimated_params_init(use_calibration);
end;

varobs P_obs Y_obs junk2;

observation_trends;
P_obs (g_p);
Y_obs (g_y);
end;

// estimated_params_init(use_calibration);
// end;

options_.plot_priors=0;
