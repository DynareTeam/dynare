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
corr gp_obs, gy_obs,0;
end;

@#define mode_file_name="fs2000_corr_ME_missing_mode"
@#define data_file_name="fsdat_simul_corr_ME_missing"

@#include "fs2000_estimation_check.inc" 