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
stderr gp_obs, inv_gamma_pdf, 0.001, inf;
end;

varobs gp_obs gy_obs;

estimation(order=1,datafile='../fsdat_simul', nobs=192, loglinear, mh_replic=2000, mh_nblocks=1, mh_jscale=0.8,forecast=8,smoother,filtered_vars,filter_step_ahead=[1:2],filter_decomposition)  m P c e W R k d y gy_obs;

if size(oo_.PointForecast.deciles.gy_obs,1)~=9
    error('Number of deciles must be 9')
end

if size(oo_.PointForecast.deciles.gy_obs,2)~=options_.forecast || ...
        size(oo_.PointForecast.Mean.gy_obs,1)~=options_.forecast || ...
        size(oo_.PointForecast.Median.gy_obs,1)~=options_.forecast || ...
        size(oo_.PointForecast.Var.gy_obs,1)~=options_.forecast || ...
        size(oo_.PointForecast.HPDinf.gy_obs,1)~=options_.forecast || ...
        size(oo_.PointForecast.HPDsup.gy_obs,1)~=options_.forecast || ...
        size(oo_.MeanForecast.deciles.gy_obs,2)~=options_.forecast || ...
        size(oo_.MeanForecast.Mean.gy_obs,1)~=options_.forecast || ...
        size(oo_.MeanForecast.Median.gy_obs,1)~=options_.forecast || ...
        size(oo_.MeanForecast.Var.gy_obs,1)~=options_.forecast || ...
        size(oo_.MeanForecast.HPDinf.gy_obs,1)~=options_.forecast || ...
        size(oo_.MeanForecast.HPDsup.gy_obs,1)~=options_.forecast 
    error('Forecasts have the wrong length')
end

if size(oo_.UpdatedVariables.deciles.gy_obs,2)~=options_.nobs || ...
        size(oo_.UpdatedVariables.Mean.gy_obs,1)~=options_.nobs || ...
        size(oo_.UpdatedVariables.Median.gy_obs,1)~=options_.nobs || ...
        size(oo_.UpdatedVariables.Var.gy_obs,1)~=options_.nobs || ...
        size(oo_.UpdatedVariables.HPDinf.gy_obs,1)~=options_.nobs || ...
        size(oo_.UpdatedVariables.HPDsup.gy_obs,1)~=options_.nobs
    error('Updated Variables have the wrong length')
end

if size(oo_.Filtered_Variables_1_step_ahead.deciles.gy_obs,2)~=options_.nobs || ...
        size(oo_.Filtered_Variables_1_step_ahead.Mean.gy_obs,1)~=options_.nobs || ...
        size(oo_.Filtered_Variables_1_step_ahead.Median.gy_obs,1)~=options_.nobs || ...
        size(oo_.Filtered_Variables_1_step_ahead.Var.gy_obs,1)~=options_.nobs || ...
        size(oo_.Filtered_Variables_1_step_ahead.HPDinf.gy_obs,1)~=options_.nobs || ...
        size(oo_.Filtered_Variables_1_step_ahead.HPDsup.gy_obs,1)~=options_.nobs
    error('Filtered_Variables_1_step_ahead have the wrong length')
end


if size(oo_.SmoothedVariables.deciles.gy_obs,2)~=options_.nobs || ...
        size(oo_.SmoothedVariables.Mean.gy_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedVariables.Median.gy_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedVariables.Var.gy_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedVariables.HPDinf.gy_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedVariables.HPDsup.gy_obs,1)~=options_.nobs
    error('Smoothed Variables have the wrong length')
end

if size(oo_.SmoothedShocks.deciles.e_a,2)~=options_.nobs || ...
        size(oo_.SmoothedShocks.Mean.e_a,1)~=options_.nobs || ...
        size(oo_.SmoothedShocks.Median.e_a,1)~=options_.nobs || ...
        size(oo_.SmoothedShocks.Var.e_a,1)~=options_.nobs || ...
        size(oo_.SmoothedShocks.HPDinf.e_a,1)~=options_.nobs || ...
        size(oo_.SmoothedShocks.HPDsup.e_a,1)~=options_.nobs
    error('Smoothed Shocks have the wrong length')
end

if size(oo_.FilteredVariables.deciles.gy_obs,2)~=options_.nobs || ...
        size(oo_.FilteredVariables.Mean.gy_obs,1)~=options_.nobs || ...
        size(oo_.FilteredVariables.Median.gy_obs,1)~=options_.nobs || ...
        size(oo_.FilteredVariables.Var.gy_obs,1)~=options_.nobs || ...
        size(oo_.FilteredVariables.HPDinf.gy_obs,1)~=options_.nobs || ...
        size(oo_.FilteredVariables.HPDsup.gy_obs,1)~=options_.nobs
    error('Filtered Variables have the wrong length')
end

if size(oo_.SmoothedMeasurementErrors.deciles.SE_EOBS_gp_obs,2)~=options_.nobs || ...
        size(oo_.SmoothedMeasurementErrors.Mean.SE_EOBS_gp_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedMeasurementErrors.Median.SE_EOBS_gp_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedMeasurementErrors.Var.SE_EOBS_gp_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedMeasurementErrors.HPDinf.SE_EOBS_gp_obs,1)~=options_.nobs || ...
        size(oo_.SmoothedMeasurementErrors.HPDsup.SE_EOBS_gp_obs,1)~=options_.nobs
    error('SmoothedMeasurementErrors have the wrong length')
end

