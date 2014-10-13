/*
 * This file replicates the estimation of the cash in advance model described
 * Frank Schorfheide (2000): "Loss function-based evaluation of DSGE models",
 * Journal of Applied Econometrics, 15(6), 645-670.
 *
 * The data are in file "fsdat_simul.m", and have been artificially generated.
 * They are therefore different from the original dataset used by Schorfheide.
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
 * Copyright (C) 2004-2010 Dynare Team
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

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
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

steady;

check;

estimated_params;
alp, beta_pdf, 0.356, 0.02;
bet, beta_pdf, 0.993, 0.002;
gam, normal_pdf, 0.0085, 0.003;
mst, normal_pdf, 1.0002, 0.007;
rho, beta_pdf, 0.129, 0.223;
psi, beta_pdf, 0.65, 0.05;
del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
end;

varobs gp_obs gy_obs;
options_.plot_priors=0;
options_.lyapunov_fp = 0;
options_.lyapunov_db = 0;
options_.lyapunov_srs = 0;

estimation(lyapunov=doubling,mode_compute=0,mode_file=fs2000_mode,order=1,datafile=fsdat_simul, nobs=192, loglinear, mh_replic=0, mh_nblocks=1, mh_jscale=0.8,nograph);

options_.lyapunov_fp = 0;
options_.lyapunov_db = 0;
options_.lyapunov_srs = 0;
estimation(lyapunov=square_root_solver,mode_compute=0,mode_file=fs2000_mode,order=1,datafile=fsdat_simul, nobs=192, loglinear, mh_replic=0, mh_nblocks=1, mh_jscale=0.8,nograph);

options_.lyapunov_fp = 0;
options_.lyapunov_db = 0;
options_.lyapunov_srs = 0;
estimation(lyapunov=fixed_point,mode_compute=0,mode_file=fs2000_mode,order=1,datafile=fsdat_simul, nobs=192, loglinear, mh_replic=0, mh_nblocks=1, mh_jscale=0.8,nograph);

//**************** Do pure Lyapunov tests***********
verbatim;
clear all
options_.lyapunov_complex_threshold = 1e-15;
options_.qz_zero_threshold = 1e-6;
options_.qz_criterium=1-options_.qz_zero_threshold;
options_.lyapunov_fixed_point_tol = 1e-10;
options_.lyapunov_doubling_tol = 1e-16;

n_small=8;
m_small=10;
T_small=randn(n_small,n_small);
T_small=0.99*T_small/max(abs(eigs(T_small)));
tmp2=randn(m_small,m_small);
Q_small=tmp2*tmp2';
R_small=randn(n_small,m_small);

n_large=9;
m_large=11;
T_large=randn(n_large,n_large);
T_large=0.99*T_large/max(abs(eigs(T_large)));
tmp2=randn(m_large,m_large);
Q_large=tmp2*tmp2';
R_large=randn(n_large,m_large);

% DynareOptions.lyapunov_fp == 1
Pstar1_small = lyapunov_symm(T_small,R_small*Q_small*R_small',options_.lyapunov_fixed_point_tol,options_.lyapunov_fixed_point_tol,3);
Pstar1_large = lyapunov_symm(T_large,R_large*Q_large*R_large',options_.lyapunov_fixed_point_tol,options_.lyapunov_fixed_point_tol,3);
% Dynareoptions.lyapunov_db == 1
Pstar2_small = disclyap_fast(T_small,R_small*Q_small*R_small',options_.lyapunov_doubling_tol);
Pstar2_large = disclyap_fast(T_large,R_large*Q_large*R_large',options_.lyapunov_doubling_tol);
% Dynareoptions.lyapunov_srs == 1
Pstar3_small = lyapunov_symm(T_small,Q_small,options_.lyapunov_fixed_point_tol,options_.lyapunov_complex_threshold,4,R_small);
Pstar3_large = lyapunov_symm(T_large,Q_large,options_.lyapunov_fixed_point_tol,options_.lyapunov_complex_threshold,4,R_large);
% Standard
Pstar4_small = lyapunov_symm(T_small,R_small*Q_small*R_small',options_.qz_criterium,options_.lyapunov_complex_threshold);
Pstar4_large = lyapunov_symm(T_large,R_large*Q_large*R_large',options_.qz_criterium,options_.lyapunov_complex_threshold);


if max(max(abs(Pstar1_small-Pstar2_small)))>1e-8
    error('Results do not match')
end
if max(max(abs(Pstar1_small-Pstar3_small)))>1e-8
    error('Results do not match')
end
if max(max(abs(Pstar1_small-Pstar4_small)))>1e-8
    error('Results do not match')
end


if max(max(abs(Pstar1_large-Pstar2_large)))>1e-8
    error('Results do not match')
end
if max(max(abs(Pstar1_large-Pstar3_large)))>1e-8
    error('Results do not match')
end
if max(max(abs(Pstar1_large-Pstar4_large)))>1e-8
    error('Results do not match')
end
    
oo_=[];
M_=[];
end;