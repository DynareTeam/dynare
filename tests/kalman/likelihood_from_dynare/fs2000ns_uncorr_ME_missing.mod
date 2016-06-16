@#include "fs2000ns_model.inc" 

estimated_params;
alp, 0.33;
gam,  0.0;
del, 0.02;
stderr e_a, 0.014;
stderr e_m, 0.005;
corr e_m, e_a, 0;
stderr P_obs, 0.05;
stderr Y_obs, 0.05;
//corr gp_obs, gy_obs,0;
end;

@#define mode_file_name="fs2000ns_uncorr_ME_missing_mode"
@#define data_file_name="fs_ns_dat_simul_uncorr_ME_missing"

@#include "fs2000ns_estimation_check.inc" 
