@#include "fs2000_model.inc" 

stoch_simul(periods=200, order=1,irf=0);
temp=oo_.endo_simul;
%add measurement error
oo_.endo_simul(strmatch('gy_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('gy_obs',M_.endo_names,'exact'),:)+0.05*randn(1,size(oo_.endo_simul,2));
oo_.endo_simul(strmatch('gp_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('gp_obs',M_.endo_names,'exact'),:)+0.05*randn(1,size(oo_.endo_simul,2));
datatomfile('fsdat_simul_uncorr_ME', char('gy_obs', 'gp_obs'));
oo_.endo_simul(strmatch('gy_obs',M_.endo_names,'exact'),[7,199])=NaN;
oo_.endo_simul(strmatch('gp_obs',M_.endo_names,'exact'),[151,199])=NaN;
datatomfile('fsdat_simul_uncorr_ME_missing', char('gy_obs', 'gp_obs'));
shock_mat=chol([1 0.5; 0.5 1])*0.05*randn(2,size(oo_.endo_simul,2));
oo_.endo_simul(strmatch('gy_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('gy_obs',M_.endo_names,'exact'),:)+shock_mat(1,:);
oo_.endo_simul(strmatch('gp_obs',M_.endo_names,'exact'),:)=oo_.endo_simul(strmatch('gp_obs',M_.endo_names,'exact'),:)+shock_mat(2,:);
datatomfile('fsdat_simul_corr_ME', char('gy_obs', 'gp_obs'));
oo_.endo_simul(strmatch('gy_obs',M_.endo_names,'exact'),[7,199])=NaN;
oo_.endo_simul(strmatch('gp_obs',M_.endo_names,'exact'),[151,199])=NaN;
datatomfile('fsdat_simul_corr_ME_missing', char('gy_obs', 'gp_obs'));

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

@#define mode_file_name="fs2000_uncorr_ME_mode"
@#define data_file_name="fsdat_simul_uncorr_ME"

@#include "fs2000_estimation_check.inc" 