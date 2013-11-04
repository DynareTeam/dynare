//uses a prior mean for bet that does not allow solving the model. Then uses the calibrated value to initialize.
var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m e_f;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
end;

initval;
k = 6;
m = mst;
P = 2.25;
c = 0.45;
e = 1;
W = 4;
R = 1.02;
d = 0.85;
n = 0.19;
l = 0.86;
y = 0.6;
gy_obs = exp(gam);
gp_obs = exp(-gam);
dA = exp(gam);
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
var e_f; stderr 0.005;
corr gy_obs, gp_obs = 0.1;
end;

steady;
//check;

estimated_params;
alp, beta_pdf, 0.356, 0.02;
bet, normal_pdf, 2, 1;
gam, normal_pdf, 0.0085, 0.003;
mst, normal_pdf, 1.0002, 0.007;
rho, normal_pdf, 0.129, 0.223;
psi, beta_pdf, 0.65, 0.05;
del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
corr e_a, e_m, normal_pdf, 0, 0.007;
corr e_a, e_f, normal_pdf, 0, 0.007;
stderr gp_obs, inv_gamma_pdf, 0.008862, inf;
corr gp_obs, gy_obs, normal_pdf, 0, 0.007;
end;

varobs gp_obs gy_obs;

estimated_params_init;
//stderr e_a, 0.014000;
//stderr e_m, 0.005000;
//stderr gp_obs, 0.1;
corr e_m, e_a, 0.1;
//corr gy_obs, gp_obs, 0.000000;
//alp, 0.330000;
//bet, 0.990000;
gam, 0.003000;
mst, 1.011000;
//rho, 0.900000;
psi, 0.787000;
del, 0.020000;
end;

options_.use_calibration_initialization=1;
options_.plot_priors=0;
estimation(order=1, datafile=fsdat_simul, nobs=192, loglinear, mh_replic=2000, mh_nblocks=1, mh_jscale=0.8,prior_trunc=0);