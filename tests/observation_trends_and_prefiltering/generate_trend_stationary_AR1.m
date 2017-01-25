function generate_trend_stationary_AR1
n_periods=10000;
rho_y=0.5;
rho_p=0.5;
g_y=0.0001;
g_p=-0.0001;
const_y=2;
const_p=2;
sigma_y=0.001;
sigma_p=0.001;

orig_params=[rho_y rho_p  g_y g_p sigma_y sigma_p]';
param_names=char('rho_y','rho_p','g_y','g_p','sigma_y','sigma_p');

save orig_params_prefilter.mat orig_params param_names

orig_params=[rho_y rho_p  g_y g_p const_y const_p sigma_y sigma_p]';
param_names=char('rho_y','rho_p','g_y','g_p','const_y','const_p','sigma_y','sigma_p');

save orig_params.mat orig_params param_names

jumping_covariance=diag([1e-8; 1e-8; 1e-16; 1e-16; 1e-8; 1e-8; 1e-12; 1e-12;])^-1;
save MCMC_jump_covar.mat jumping_covariance
jumping_covariance=diag([1e-8; 1e-8; 1e-16; 1e-16; 1e-12; 1e-12;])^-1;
save MCMC_jump_covar_prefilter.mat jumping_covariance

%% data without constant
log_P=zeros(1,n_periods);
log_Y=zeros(1,n_periods);
junk2_orig=zeros(1,n_periods);
for ii=2:n_periods
    log_P(ii)=rho_p*log_P(ii-1)+sigma_p*randn;
    log_Y(ii)=rho_y*log_Y(ii-1)+sigma_y*randn;
    junk2_orig(ii)=0.9*junk2_orig(ii-1)+randn;
end
%add trend
log_P=log_P+g_p*(1:n_periods);
log_Y=log_Y+g_y*(1:n_periods);

Y_obs=exp(log_Y);
P_obs=exp(log_P);
junk2=exp(junk2_orig);
save Exp_AR1_trend_data_no_constant.mat Y_obs P_obs junk2
% 
% [b_p,~,~,~,stats_p] = regress(log(P_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' log(P_obs(1:end-1)')]);
% [b_y,~,~,~,stats_y] = regress(log(Y_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' log(Y_obs(1:end-1)')]);

Y_obs=log_Y;
P_obs=log_P;
junk2=junk2_orig;  
save AR1_trend_data_no_constant.mat Y_obs P_obs junk2

% [b_p,~,~,~,stats_p] = regress((P_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' (P_obs(1:end-1)')]);
% [b_y,~,~,~,stats_y] = regress((Y_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' (Y_obs(1:end-1)')]);

%% data with constant
log_P=zeros(1,n_periods);
log_Y=zeros(1,n_periods);
log_P(1,1)=const_p;
log_Y(1,1)=const_y;
for ii=2:n_periods
    log_P(ii)=(1-rho_p)*const_p+rho_p*log_P(ii-1)+sigma_p*randn;
    log_Y(ii)=(1-rho_y)*const_y+rho_y*log_Y(ii-1)+sigma_y*randn;
end
%add trend
log_P=log_P+g_p*(1:n_periods);
log_Y=log_Y+g_y*(1:n_periods);

Y_obs=exp(log_Y);
P_obs=exp(log_P);
junk2=exp(junk2_orig);
save Exp_AR1_trend_data_with_constant.mat Y_obs P_obs junk2

% [b,bint,r,rint,stats] = regress(log(P_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' log(P_obs(1:end-1)')]);
% [b,bint,r,rint,stats] = regress(log(Y_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' log(Y_obs(1:end-1)')]);

Y_obs=log_Y;
P_obs=log_P;
junk2=junk2_orig;  
save AR1_trend_data_with_constant.mat Y_obs P_obs junk2

% [b_p,~,~,~,stats_p] = regress((P_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' (P_obs(1:end-1)')]);
% [b_y,~,~,~,stats_y] = regress((Y_obs(2:end))',[ones(n_periods-1,1) (2:n_periods)' (Y_obs(1:end-1)')]);


