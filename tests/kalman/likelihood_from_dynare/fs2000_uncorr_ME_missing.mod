/*
 * This file is based on the cash in advance model described
 * Frank Schorfheide (2000): "Loss function-based evaluation of DSGE models",
 * Journal of Applied Econometrics, 15(6), 645-670.
 *
 * The equations are taken from J. Nason and T. Cogley (1994): "Testing the
 * implications of long-run neutrality for monetary business cycle models",
 * Journal of Applied Econometrics, 9, S37-S70.
 * Note that there is an initial minus sign missing in equation (A1), p. S63.
 *
 * This implementation was written by Michel Juillard. Please note that the
 * following copyright notice only applies to this Dynare implementation of the
 * model.
 */

/*
 * Copyright (C) 2004-2013 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del theta;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;
theta=0;

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

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
end;

varobs gp_obs gy_obs;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

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

options_.TeX=1;
options_.debug=1;

%%default
estimation(kalman_algo=0,mode_compute=4,order=1,datafile=fsdat_simul_uncorr_ME_missing,smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_0=oo_.likelihood_at_initial_parameters;
%%Multivariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=1,mode_file=fs2000_uncorr_ME_missing_mode,mode_compute=0,order=1,datafile=fsdat_simul_uncorr_ME_missing,smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_1=oo_.likelihood_at_initial_parameters;
%%Univariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=3,mode_file=fs2000_uncorr_ME_missing_mode,mode_compute=0,order=1,datafile=fsdat_simul_uncorr_ME_missing,smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_3=oo_.likelihood_at_initial_parameters;
%%Diffuse Multivariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=2,mode_file=fs2000_uncorr_ME_missing_mode,mode_compute=0,order=1,datafile=fsdat_simul_uncorr_ME_missing,smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_2=oo_.likelihood_at_initial_parameters;
%%Diffuse univariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=4,mode_file=fs2000_uncorr_ME_missing_mode,mode_compute=0,order=1,datafile=fsdat_simul_uncorr_ME_missing,smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_4=oo_.likelihood_at_initial_parameters;

if max(abs([fval_algo_2,fval_algo_3,fval_algo_4]-fval_algo_1))>1e-6
    error('Likelihoods do not match')
end