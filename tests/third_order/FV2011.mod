/*
 * This file replicates the IRFs of the small open economy model described in
 * Jesús Fernández-Villaverde, Pablo Guerrón-Quintana,
 * Juan F. Rubio-Ramírez, and Martin Uribe (2011): "Risk Matters",
 * American Economic Review 101 (October 2011): 2530–2561.
 *
 * This implementation was written by Benjamin Born and Johannes Pfeifer. Please 
 * note that the following copyright notice only applies to this Dynare 
 * implementation of the
 * model.
 */

/*
 * Copyright (C) 2013 Dynare Team
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

var sigma_r sigma_tb eps_r eps_tb X D K lambda C H Y I phi r;
varexo u_sigma_r u_sigma_tb u_r u_tb u_x ;
predetermined_variables K D;

parameters r_bar rho_eps_r sigma_r_bar rho_sigma_r eta_r 
           rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb
           delta alppha nu rho_x betta
           Phi phipar sigma_x D_bar omega eta;

// Calibration from Table 3
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94; 
eta_r=0.46;

// Calibration from Table 4
rho_eps_tb=0.95;
sigma_tb_bar=-8.06; %8.05 in paper, but 8.06 in code
rho_sigma_tb=0.94;
eta_tb=0.13;
nu=5; %inverse of the elasticity of intertemporal substitution
eta=1000; %elasticity of labor to wages
delta=0.014; %depreciation
alppha = 0.32; % capital income share
rho_x=0.95; %autocorrelation TFP


// Calibration from Table 6 //
r_bar=log(0.02);
betta = 1/(1+exp(r_bar)); %discount factor
Phi=0.001; %debt elasticity
D_bar= 4; % steady state debt
phipar=95; % capital adjustment costs
sigma_x=log(0.015); 

omega=1;

model;
(exp(C))^-nu=exp(lambda);
exp(lambda)/(1+exp(r))=exp(lambda)*Phi*((D(+1))-D_bar)+betta*exp(lambda(+1));
-exp(phi)+betta*((1-delta)*exp(phi(+1))+alppha*exp(Y(+1))/exp(K(+1))*exp(lambda(+1)))=0;
omega*exp(H)^eta=(1-alppha)*exp(Y)/exp(H)*exp(lambda);
exp(phi)*(1-phipar/2*((exp(I)-exp(I(-1)))/exp(I(-1)))^2-phipar*exp(I)/exp(I(-1))*((exp(I)-exp(I(-1)))/exp(I(-1))))+betta*exp(phi(+1))*phipar*(exp(I(+1))/exp(I))^2*((exp(I(+1))-exp(I))/exp(I))=exp(lambda);

exp(Y)=exp(K)^alppha*(exp(X)*exp(H))^(1-alppha);
X=rho_x*X(-1)+exp(sigma_x)*u_x;
exp(K(+1))=(1-delta)*exp(K)+(1-phipar/2*(exp(I)/exp(I(-1))-1)^2)*exp(I);
exp(Y)-exp(C)-exp(I)=(D)-(D(+1))/(1+exp(r))+Phi/2*((D(+1))-D_bar)^2;
exp(r)=exp(r_bar)+eps_tb+eps_r;
eps_tb=rho_eps_tb*eps_tb(-1)+exp(sigma_tb)*u_tb;   
sigma_tb=(1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb(-1)+eta_tb*u_sigma_tb;
eps_r=rho_eps_r*eps_r(-1)+exp(sigma_r)*u_r;   
sigma_r=(1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r(-1)+eta_r*u_sigma_r;
end;

initval;
sigma_tb=sigma_tb_bar;
sigma_r=sigma_r_bar;
eps_r=0;
eps_tb=0;
D=D_bar;

% steady states taken from Mathematica code
C=0.8779486025329908;
K=3.293280327636415;
lambda=-4.389743012664954;
H=-0.0037203652717462993;
phi=-4.389743012664954;
I=-0.9754176217303792;
Y=1.0513198564588924;
r=r_bar; 
end;

shocks;
var u_x; stderr 1;
var u_r; stderr 1;
var u_tb; stderr 1;
var u_sigma_tb; stderr 1;
var u_sigma_r; stderr 1;
end;

resid(1);

options_.solve_tolf=1E-12;
steady(solve_algo=3);

check;
stoch_simul(order=3,pruning,irf=0,nocorr,nofunctions,nomoments) C I Y H r D K lambda phi; 

comparison_policy_functions_dynare_mathematica;