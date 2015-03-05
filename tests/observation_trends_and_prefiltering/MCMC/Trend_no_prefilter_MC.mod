var Y_obs P_obs junk1 junk2;
varexo e_y e_p;

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
Y_obs = (1-rho_y)*const_y + rho_y*Y_obs(-1)+sigma_y*e_y;
P_obs = (1-rho_p)*const_p + rho_p*P_obs(-1)+sigma_p*e_p;
junk1 = 0.9*junk1(+1);
junk2 = 0.9*junk2(-1);
end;

steady_state_model;
Y_obs = const_y;
P_obs = const_p;
junk1=0;
junk2=0;
end;

shocks;
var e_p; stderr 1;
var e_y; stderr 1;
end;

steady;
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

estimated_params_init(use_calibration);
end;

options_.plot_priors=0;
estimation(order=1,datafile='../AR1_trend_data_with_constant',mh_replic=2000,
            mode_compute=4,first_obs=1,smoother,mh_nblocks=1,mh_jscale=0.3,
            mcmc_jumping_covariance='MCMC_jump_covar',forecast=100,prefilter=0) P_obs Y_obs junk2;
            
load('../AR1_trend_data_with_constant');
loaded_par=load('../orig_params');
if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end
y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',loaded_par.param_names,'exact'));
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',loaded_par.param_names,'exact'));

if abs(oo_.PointForecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.PointForecast.Mean.P_obs(end)- p_forecast_100_periods)>5e-4
    error('Mean Point Forecasts do not match')
end
if abs(oo_.PointForecast.Median.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.PointForecast.Median.P_obs(end)- p_forecast_100_periods)>5e-4
    error('Median Point Forecasts do not match')
end

if abs(oo_.MeanForecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.MeanForecast.Mean.P_obs(end)- p_forecast_100_periods)>2e-4
    error('Mean Mean Forecasts do not match')
end
if abs(oo_.MeanForecast.Median.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.MeanForecast.Median.P_obs(end)- p_forecast_100_periods)>2e-3
    error('Median Mean Forecasts do not match')
end

if abs(mean(oo_.SmoothedShocks.Mean.e_y))>1e-2 || abs(mean(oo_.SmoothedShocks.Mean.e_p))>1e-2 || abs(mean(oo_.SmoothedShocks.Median.e_y))>1e-2 || abs(mean(oo_.SmoothedShocks.Median.e_p))>1e-2
    error('Residuals are not mean 0')
end