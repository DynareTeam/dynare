var Y_obs P_obs junk1 junk2;
varexo e_y e_p;

parameters rho_y rho_p  g_y g_p sigma_y sigma_p;

rho_y=0.5;
rho_p=0.5;
g_y=0.0001;
g_p=-0.0001;
sigma_y=0.001;
sigma_p=0.001;


model;
Y_obs = rho_y*Y_obs(-1)+sigma_y*e_y;
P_obs = rho_p*P_obs(-1)+sigma_p*e_p;
junk1 = 0.9*junk1(+1);
junk2 = 0.9*junk2(-1);

end;

steady_state_model;
Y_obs = 0;
P_obs = 0;
junk1=0;
junk2=0;
end;

shocks;
var e_p; stderr 1;
var e_y; stderr 1;
end;

steady(nocheck);
check;

estimated_params;
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

estimated_params_init(use_calibration);
end;

options_.plot_priors=0;

estimation(order=1,datafile='../AR1_trend_data_with_constant',mh_replic=0,mode_compute=4,
first_obs=1,diffuse_filter,smoother,forecast=100,prefilter=1) P_obs Y_obs;

load('../AR1_trend_data_with_constant');
loaded_par=load('../orig_params');

if max(abs((M_.params-loaded_par.orig_params([1:4,7:8]))./loaded_par.orig_params([1:4,7:8])))>0.03
    error('Parameter estimates do not match')
end

y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',loaded_par.param_names,'exact'))
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',loaded_par.param_names,'exact'))

if abs(oo_.forecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.forecast.Mean.P_obs(end)- p_forecast_100_periods)>1e-4
    error('Forecasts do not match')
end

if abs(mean(oo_.SmoothedShocks.e_y))>1e-4 || abs(mean(oo_.SmoothedShocks.e_p))>1e-4
    error('Residuals are not mean 0')
end