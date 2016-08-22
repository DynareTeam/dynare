@#include "fs2000ns_model.inc" 

stoch_simul(periods=200, order=1,irf=0);
temp=oo_.endo_simul;
%add measurement error
oo_.endo_simul(strmatch('Y_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('Y_obs',M_.endo_names,'exact'),:)+0.05*randn(1,size(oo_.endo_simul,2));
oo_.endo_simul(strmatch('P_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('P_obs',M_.endo_names,'exact'),:)+0.05*randn(1,size(oo_.endo_simul,2));
datatomfile('fs_ns_dat_simul_uncorr_ME', char('Y_obs', 'P_obs'));
oo_.endo_simul(strmatch('Y_obs',M_.endo_names,'exact'),[7,199])=NaN;
oo_.endo_simul(strmatch('P_obs',M_.endo_names,'exact'),[151,199])=NaN;
datatomfile('fs_ns_dat_simul_uncorr_ME_missing', char('Y_obs', 'P_obs'));
shock_mat=chol([1 0.5; 0.5 1])*0.05*randn(2,size(oo_.endo_simul,2));
oo_.endo_simul(strmatch('Y_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('Y_obs',M_.endo_names,'exact'),:)+shock_mat(1,:);
oo_.endo_simul(strmatch('P_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('P_obs',M_.endo_names,'exact'),:)+shock_mat(2,:);
datatomfile('fs_ns_dat_simul_corr_ME', char('Y_obs', 'P_obs'));
oo_.endo_simul(strmatch('Y_obs',M_.endo_names,'exact'),[7,199])=NaN;
oo_.endo_simul(strmatch('P_obs',M_.endo_names,'exact'),[151,199])=NaN;
datatomfile('fs_ns_dat_simul_corr_ME_missing', char('Y_obs', 'P_obs'));

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

@#define mode_file_name="fs2000ns_uncorr_ME_mode"
@#define data_file_name="fs_ns_dat_simul_uncorr_ME"

@#include "fs2000ns_estimation_check.inc" 

