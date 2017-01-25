       
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

if isoctave
    tol_crit=3e-2 ;
else
    tol_crit=2e-2 ;    
end

if abs(corr(oo_.FilteredVariables.Mean.Y_obs(2:end-1)-Y_obs(3:end)',oo_.FilteredVariables.Mean.Y_obs(1:end-2)-Y_obs(2:end-1)'))>2e-2 ||...
    abs(corr(oo_.FilteredVariables.Mean.P_obs(2:end-1)-P_obs(3:end)',oo_.FilteredVariables.Mean.P_obs(1:end-2)-P_obs(2:end-1)'))>2e-2 ||...
    abs(corr(oo_.FilteredVariables.Mean.junk2(2:end-1)-junk2(3:end)',oo_.FilteredVariables.Mean.junk2(1:end-2)-junk2(2:end-1)'))>tol_crit
    error('Filtered Variables are wrong')
end
   
if max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,1,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.P_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,2,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.Y_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,3,2:end-(options_.nk-1)))-oo_.FilteredVariables.Mean.junk2))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(2,1,3:end-options_.nk))-oo_.FilteredVariables.Mean.P_obs(3:end)))>1e-2 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(2,2,3:end-options_.nk))-oo_.FilteredVariables.Mean.Y_obs(3:end)))>1e-2 ||...
    mean(squeeze(oo_.FilteredVariablesKStepAhead(2,3,3:end-options_.nk)))>1e-1 
    error('FilteredVariablesKStepAhead is wrong')
end

if isoctave
  tol_crit=7e-4;
else
  tol_crit=5e-4;
end
        
if abs(oo_.PointForecast.Mean.Y_obs(end)- y_forecast_100_periods)>tol_crit || abs(oo_.PointForecast.Mean.P_obs(end)- p_forecast_100_periods)>tol_crit
    error('Mean Point Forecasts do not match')
end
if abs(oo_.PointForecast.Median.Y_obs(end)- y_forecast_100_periods)>tol_crit || abs(oo_.PointForecast.Median.P_obs(end)- p_forecast_100_periods)>tol_crit
    error('Median Point Forecasts do not match')
end

if abs(oo_.MeanForecast.Mean.Y_obs(end)- y_forecast_100_periods)>tol_crit || abs(oo_.MeanForecast.Mean.P_obs(end)- p_forecast_100_periods)>tol_crit
    error('Mean Mean Forecasts do not match')
end
if abs(oo_.MeanForecast.Median.Y_obs(end)- y_forecast_100_periods)>tol_crit || abs(oo_.MeanForecast.Median.P_obs(end)- p_forecast_100_periods)>5e-3
    error('Median Mean Forecasts do not match')
end

if abs(mean(oo_.SmoothedShocks.Mean.e_y))>1e-2 || abs(mean(oo_.SmoothedShocks.Mean.e_p))>1e-2 || abs(mean(oo_.SmoothedShocks.Median.e_y))>1e-2 || abs(mean(oo_.SmoothedShocks.Median.e_p))>1e-2
    error('Residuals are not mean 0')
end