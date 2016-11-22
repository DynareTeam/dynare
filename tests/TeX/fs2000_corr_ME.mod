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
 * Copyright (C) 2004-2016 Dynare Team
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

var m ${m}$ (long_name='Money Stock')
    P ${P}$ (long_name='Price Level')
    c ${c}$ (long_name='Consumption')
    e ${e}$ (long_name='exp(Tech. Shock)')
    W ${W}$ (long_name='Nominal Wage')
    R ${R}$ (long_name='Nominal Rental Rate of Capital')
    k ${k}$ (long_name='Capital')
    d ${d}$ (long_name='Deposits')
    n ${n}$ (long_name='Hours worked')
    l ${l}$ (long_name='Loans')
    gy_obs ${\Delta y^{obs}}$ (long_name='Observed growth rate of output')
    gp_obs ${\Delta m^{obs}}$ (long_name='Observed growth rate of prices')
    y  ${y}$ (long_name='Output')
    dA ${\Delta A}$ (long_name='Labor Augm. Techn. Growth Rate')
            ;
varexo e_a ${\varepsilon_a}$ (long_name='Technology shock')
    e_m ${\varepsilon_m}$ (long_name='Observed money growth rate')
        ;

parameters alp ${\alpha}$ (long_name='capital share')
        bet ${\beta}$ (long_name='discount factor')
        gam ${\gamma}$ (long_name='Average technology growth')
        mst ${\bar m}$ (long_name='Average money stock')
        rho ${\rho}$ (long_name='Autocorrelation money process')
        psi ${\psi}$ (long_name='Leisure weight in utility')
        del ${\delta}$ (long_name='depreciation')
        ;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

options_.TeX=1;

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
corr gy_obs,gp_obs = 0.5;
end;

steady;

stoch_simul(order=1,irf=20,graph_format=eps,periods=1000,contemporaneous_correlation,conditional_variance_decomposition=[1,3]);
options_.rplottype=0;
rplot y e_a gy_obs;
options_.rplottype=1;
rplot l e_m gp_obs;
options_.rplottype=2;
rplot n e_a e_m m;
stoch_simul(order=1,irf=20,graph_format=eps,periods=0,contemporaneous_correlation,conditional_variance_decomposition=[1,3]);

write_latex_original_model;
write_latex_static_model;
write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

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

estimation(order=1,datafile='../fs2000/fsdat_simul',mode_check,smoother,filter_covariance,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20,contemporaneous_correlation) m P c e W R k d y gy_obs;



estimated_params;
//alp, beta_pdf, 0.356, 0.02;
gam, normal_pdf, 0.0085, 0.003;
//del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
corr e_m, e_a, normal_pdf, 0, 0.2;
stderr gp_obs, inv_gamma_pdf, 0.001, inf;
stderr gy_obs, inv_gamma_pdf, 0.001, inf;
corr gp_obs, gy_obs,normal_pdf, 0, 0.2;
end;

write_latex_prior_table;

estimation(mode_compute=8,order=1,datafile='../fs2000/fsdat_simul',mode_check,smoother,filter_decomposition,mh_replic=4000, mh_nblocks=1, mh_jscale=0.8,forecast = 8,bayesian_irf,filtered_vars,filter_step_ahead=[1,3],irf=20,
        moments_varendo,contemporaneous_correlation,conditional_variance_decomposition=[1 2 4],smoothed_state_uncertainty,raftery_lewis_diagnostics) m P c e W R k d y gy_obs;

trace_plot(options_,M_,estim_params_,'PosteriorDensity',1);
trace_plot(options_,M_,estim_params_,'StructuralShock',1,'e_a')

shock_decomposition y W R;

collect_latex_files;

//identification(advanced=1,max_dim_cova_group=3,prior_mc=250);
if system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])
    error('TeX-File did not compile.')
end
