@#include "../Trend_exp_model_calib_no_prefilter_common.inc"
options_.filter_decomposition=1;

calib_smoother(datafile='../Exp_AR1_trend_data_with_constant',prefilter=0,loglinear,
//         filter_decomposition,
        filtered_vars, filter_step_ahead = [1,2,4]) P_obs Y_obs junk2;
load('../Exp_AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 


loaded_par=load('../orig_params');
if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameters do not match')
end

@#include "../Trend_diagnostics_calib_common.inc"