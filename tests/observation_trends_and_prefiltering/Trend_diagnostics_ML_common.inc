if abs(oo_.forecast.Mean.Y_obs(end)- y_forecast_100_periods)>1e-4 || abs(oo_.forecast.Mean.P_obs(end)- p_forecast_100_periods)>1e-4
    error('Forecasts do not match')
end

if max(abs(oo_.SmoothedVariables.Y_obs-Y_obs'))>1e-5 ||...
    max(abs(oo_.SmoothedVariables.P_obs-P_obs'))>1e-5 || ...
    max(abs(oo_.SmoothedVariables.junk2-junk2'))>1e-5
    error('Smoothed Variables are wrong')
end

if max(abs(oo_.UpdatedVariables.Y_obs-Y_obs'))>1e-5 ||...
    max(abs(oo_.UpdatedVariables.P_obs-P_obs'))>1e-5 || ...
    max(abs(oo_.UpdatedVariables.junk2-junk2'))>1e-5
    error('Updated Variables are wrong')
end

if mean(abs(oo_.FilteredVariables.Y_obs(1:end-1)-Y_obs(2:end)'))>1e-3 ||...
    mean(abs(oo_.FilteredVariables.P_obs(1:end-1)-P_obs(2:end)'))>1e-3 
    error('Smoothed Variables are wrong')
end

if ~isoctave
  if abs(corr(oo_.FilteredVariables.Y_obs(2:end-1)-Y_obs(3:end)',oo_.FilteredVariables.Y_obs(1:end-2)-Y_obs(2:end-1)'))>2e-2 ||...
      abs(corr(oo_.FilteredVariables.P_obs(2:end-1)-P_obs(3:end)',oo_.FilteredVariables.P_obs(1:end-2)-P_obs(2:end-1)'))>2e-1 ||...
      abs(corr(oo_.FilteredVariables.junk2(2:end-1)-junk2(3:end)',oo_.FilteredVariables.junk2(1:end-2)-junk2(2:end-1)'))>2e-2 
      error('Filtered Variables are wrong')
  end
else
  if abs(corr(oo_.FilteredVariables.Y_obs(2:end-1)-Y_obs(3:end)',oo_.FilteredVariables.Y_obs(1:end-2)-Y_obs(2:end-1)'))>2e-2 ||...
      abs(corr(oo_.FilteredVariables.P_obs(2:end-1)-P_obs(3:end)',oo_.FilteredVariables.P_obs(1:end-2)-P_obs(2:end-1)'))>2e-1 ||...
      abs(corr(oo_.FilteredVariables.junk2(2:end-1)-junk2(3:end)',oo_.FilteredVariables.junk2(1:end-2)-junk2(2:end-1)'))>3e-2 
      error('Filtered Variables are wrong')
  end
end
   
Y_pos=strmatch('Y_obs',M_.endo_names,'exact');
P_pos=strmatch('P_obs',M_.endo_names,'exact');
junk2_pos=strmatch('junk2',M_.endo_names,'exact');

if max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,Y_pos,2:end-(options_.nk-1)))-oo_.FilteredVariables.Y_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,P_pos,2:end-(options_.nk-1)))-oo_.FilteredVariables.P_obs))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(1,junk2_pos,2:end-(options_.nk-1)))-oo_.FilteredVariables.junk2))>1e-5 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(2,Y_pos,3:end-options_.nk))-oo_.FilteredVariables.Y_obs(3:end)))>1e-2 ||...
    max(abs(squeeze(oo_.FilteredVariablesKStepAhead(2,P_pos,3:end-options_.nk))-oo_.FilteredVariables.P_obs(3:end)))>1e-2 ||...
    mean(squeeze(oo_.FilteredVariablesKStepAhead(2,3,3:end-options_.nk)))>1e-1 ||...
    mean(squeeze(oo_.FilteredVariablesKStepAhead(2,4,3:end-options_.nk)))>1e-1
    error('FilteredVariablesKStepAhead is wrong')
end
   
if abs(mean(oo_.SmoothedShocks.e_y))>1e-1 || abs(mean(oo_.SmoothedShocks.e_p))>1e-1
    error('Residuals are not mean 0')
end