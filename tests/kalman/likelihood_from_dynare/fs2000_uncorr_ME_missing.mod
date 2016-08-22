@#include "fs2000_model.inc" 

estimated_params;
alp, 0.356;
gam,  0.0085;
del, 0.01;
stderr e_a, 0.035449;
stderr e_m, 0.008862;
corr e_m, e_a, 0;
stderr gp_obs, 1;
stderr gy_obs, 1;
//corr gp_obs, gy_obs,0;
end;

@#define mode_file_name="fs2000_uncorr_ME_missing_mode"
@#define data_file_name="fsdat_simul_uncorr_ME_missing"

@#include "fs2000_estimation_check.inc" 

%%Multivariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=1,fast_kalman_filter,mode_file=@{mode_file_name},mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_5=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,1)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,1)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,1)=cell2mat(struct2cell(oo_.SmoothedVariables));