var Y_obs P_obs junk1 junk2;
varexo e_y e_p;

parameters rho_y rho_p  g_y g_p const_y const_p sigma_y sigma_p;

rho_y=0.5;
rho_p=0.5;
g_y=0.001;
g_p=-0.001;
const_y=2;
const_p=2;
sigma_y=0.001;
sigma_p=0.001;

model;
Y_obs = exp(const_y)^(1-rho_y)*Y_obs(-1)^rho_y*exp(sigma_y*e_y);
P_obs = exp(const_p)^(1-rho_p)*P_obs(-1)^rho_p*exp(sigma_p*e_p);
junk1 = (junk1(+1))^0.9;
junk2 = (junk2(-1))^0.9;
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

varobs P_obs Y_obs;

observation_trends;
P_obs (g_p);
Y_obs (g_y);
end;

options_.plot_priors=0;

estimation(order=1,datafile='../Exp_AR1_trend_data_with_constant',mh_replic=0,
    mode_file=Trend_loglinear_no_prefilter_mode,
    mode_compute=4,first_obs=1,loglinear,smoother,forecast=100,prefilter=0) P_obs Y_obs;
load('../Exp_AR1_trend_data_with_constant');
loaded_par=load('../orig_params');
if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end

y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',M_.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',M_.param_names,'exact'));
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',M_.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',M_.param_names,'exact'));

if abs(oo_.forecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.forecast.Mean.P_obs(end)- p_forecast_100_periods)>1e-4
    error('Forecasts do not match')
end

if abs(mean(oo_.SmoothedShocks.e_y))>1e-1 || abs(mean(oo_.SmoothedShocks.e_p))>1e-1
    error('Residuals are not mean 0')
end