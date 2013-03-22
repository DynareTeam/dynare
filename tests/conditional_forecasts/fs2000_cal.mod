// See fs2000.mod in the examples/ directory for details on the model

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

//model(block, bytecode);
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
end;

steady;

check;

//stoch_simul(irf=0);

conditional_forecast_paths;
var gp_obs;
periods 1 2:5;
//values  0.05;
//values  0.98 1.00797;
values  0.98 0.99;
//expectation perfect_foresight;
var gy_obs;
periods  1  2  3:5;
//values   0.01 -0.02 0;
//values   0.85 0.85 0.95;
values   0.95 0.95 0.99;
//expectation perfect_foresight;
end;

options_.stack_solve_algo = 0;
options_.maxit_ = 50;

conditional_forecast(parameter_set=calibration, controlled_varexo=(e_m,e_a), simulation_type = deterministic);

/*shocks;
var e_a;
periods 1 2 3 4 5;
values -0.0109 -0.0122 -0.0137 -0.0154 -0.0173;
var e_m;
periods 1 2 3 4 5;
values -0.1242 -0.0386 -0.0392 -0.0398 -0.0405;
end;
simul(periods=40);*/
rplot gy_obs;
rplot gp_obs;
//if ~(exist('OCTAVE_VERSION') && octave_ver_less_than('3.4.0'))
//plot_conditional_forecast(periods=10) gy_obs gp_obs;
//end
