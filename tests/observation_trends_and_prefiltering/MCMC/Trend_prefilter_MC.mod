var Y_obs P_obs junk1 junk2;
varexo e_y e_p eps_junk;

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
junk2 = 0.9*junk2(-1)+eps_junk;

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
var eps_junk; stderr 1;
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

estimation(order=1,datafile='../AR1_trend_data_with_constant',mh_replic=2000,mode_compute=4,
first_obs=1,smoother,prefilter=1,
mh_nblocks=1,mh_jscale=1e-4,
filtered_vars, filter_step_ahead = [1,2,4],
mcmc_jumping_covariance='MCMC_jump_covar_prefilter',forecast=100) P_obs Y_obs junk2;

load('../AR1_trend_data_with_constant');
loaded_par=load('../orig_params_prefilter');

if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end
loaded_par_full=load('../orig_params');
y_forecast_100_periods=loaded_par_full.orig_params(strmatch('const_y',loaded_par_full.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par_full.orig_params(strmatch('g_y',loaded_par_full.param_names,'exact'));
p_forecast_100_periods=loaded_par_full.orig_params(strmatch('const_p',loaded_par_full.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par_full.orig_params(strmatch('g_p',loaded_par_full.param_names,'exact'));


if max(abs(oo_.SmoothedVariables.Mean.Y_obs-Y_obs'))>1e-5 ||...
    max(abs(oo_.SmoothedVariables.Mean.P_obs-P_obs'))>1e-5 || ...
    max(abs(oo_.SmoothedVariables.Mean.junk2-junk2'))>1e-5
    error('Smoothed Variables are wrong')
end

if max(abs(oo_.UpdatedVariables.Mean.Y_obs-Y_obs'))>1e-5 ||...
    max(abs(oo_.UpdatedVariables.Mean.P_obs-P_obs'))>1e-5 || ...
    max(abs(oo_.UpdatedVariables.Mean.junk2-junk2'))>1e-5
    error('Updated Variables are wrong')
end

if mean(abs(oo_.FilteredVariables.Mean.Y_obs(1:end-1)-Y_obs(2:end)'))>1e-3 ||...
    mean(abs(oo_.FilteredVariables.Mean.P_obs(1:end-1)-P_obs(2:end)'))>1e-3 
    error('Filtered Variables are wrong')
end

if abs(corr(oo_.FilteredVariables.Mean.Y_obs(2:end-1)-Y_obs(3:end)',oo_.FilteredVariables.Mean.Y_obs(1:end-2)-Y_obs(2:end-1)'))>2e-2 ||...
    abs(corr(oo_.FilteredVariables.Mean.P_obs(2:end-1)-P_obs(3:end)',oo_.FilteredVariables.Mean.P_obs(1:end-2)-P_obs(2:end-1)'))>2e-2 ||...
    abs(corr(oo_.FilteredVariables.Mean.junk2(2:end-1)-junk2(3:end)',oo_.FilteredVariables.Mean.junk2(1:end-2)-junk2(2:end-1)'))>2e-2 
    error('Filtered Variables are wrong')
end

if max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,1,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.P_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,2,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.Y_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,3,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.junk2))>1e-5 
    error('FilteredVariablesKStepAhead is wrong')
end
   
if max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,1,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.P_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,2,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.Y_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,3,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.junk2))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(2,1,3:end-options_.nk))-oo_.FilteredVariables.Mean.P_obs(3:end)))>1e-2 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(2,2,3:end-options_.nk))-oo_.FilteredVariables.Mean.Y_obs(3:end)))>1e-2 ||...
    mean(squeeze(oo_.FilteredVariablesKStepAhead(2,3,3:end-options_.nk)))>1e-1 
    error('FilteredVariablesKStepAhead is wrong')
end

            
if abs(oo_.PointForecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.PointForecast.Mean.P_obs(end)- p_forecast_100_periods)>5e-4
    error('Mean Point Forecasts do not match')
end
if abs(oo_.PointForecast.Median.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.PointForecast.Median.P_obs(end)- p_forecast_100_periods)>5e-3
    error('Median Point Forecasts do not match')
end

if abs(oo_.MeanForecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.MeanForecast.Mean.P_obs(end)- p_forecast_100_periods)>1e-4
    error('Mean Mean Forecasts do not match')
end
if abs(oo_.MeanForecast.Median.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.MeanForecast.Median.P_obs(end)- p_forecast_100_periods)>1e-4
    error('Median Mean Forecasts do not match')
end

if abs(mean(oo_.SmoothedShocks.Mean.e_y))>1e-4 || abs(mean(oo_.SmoothedShocks.Mean.e_p))>1e-4 || abs(mean(oo_.SmoothedShocks.Median.e_y))>1e-4 || abs(mean(oo_.SmoothedShocks.Median.e_p))>1e-4
    error('Residuals are not mean 0')
end