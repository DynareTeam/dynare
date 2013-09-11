/*
 * This file replicates the model studied in:
 * Aguiar, Mark and Gopinath, Gita (2004): "Emerging Market Business Cycles:
 * The Cycle is the Trend" (NBER WP 10734). It is different from version published 
 * in the Journal of Political Economy.
 * 
 * This model file is intended to show the capabilities of the Dynare macro 
 * languange. It is not intended to provide a full replication of the original 
 * paper due to some differences in model calibration. In 
 * particular, this mod-file does not calibrate the share of debt to GDP 
 * to 0.1 as this would require the use of a steady state file. Rather, the
 * absolute value of debt is set to 0.1. Given that output is close to 1 in
 * the benchmark specification, this results in only a small difference to 
 * the working paper. 
 * The mod-file reproduces Figure 4 of the working paper, which displays the
 * model response to 1 percent shock to trend and cyclical TFP.
 *
 * This implementation was written by Sébastien Villemot and Johannes Pfeifer. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2012-13 Dynare Team
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


// Set the following variable to 0 to get Cobb-Douglas utility
@#define ghh = 1
// Set the following variable to 0 to get the calibration for Canada
@#define mexico = 1

var c k y b q g l u z uc ul f c_y tb_y i_y;

varexo eps_z eps_g;

parameters mu_g sigma rho_g sigma_g delta phi psi b_star alpha rho_z sigma_z r_star beta;

// Benchmark parameter values (table 3)
@#if ghh == 1
parameters tau nu;
tau = 1.4;
nu = 1.6;
@#else
parameters gamma;
gamma = 0.36;
@#endif

alpha = 0.68;
sigma = 2;
delta = 0.03;
beta = 0.98;
psi = 0.001;
b_star = 0.1; //taken here as the steady state value of debt; in the original paper, this is the share of debt to GDP

// Estimated parameters (table 4)
@#if mexico == 1
@# if ghh == 1
mu_g = 1.006;
sigma_z = 0.0041;
rho_z = 0.94;
sigma_g = 0.0109;
rho_g = 0.72;
phi = 3.79;
@# else
mu_g = 1.005;
sigma_z = 0.0046;
rho_z = 0.94;
sigma_g = 0.025;
rho_g = 0.06;
phi = 2.82;
@# endif
@#else
// Canada
@# if ghh == 1
mu_g = 1.007;
sigma_z = 0.0057;
rho_z = 0.88;
sigma_g = 0.0014;
rho_g = 0.94;
phi = 2.63;
@# else
mu_g = 1.007;
sigma_z = 0.0072;
rho_z = 0.96;
sigma_g = 0.0044;
rho_g = 0.50;
phi = 3.76;
@# endif
@#endif

@#if ghh == 1
r_star = mu_g^sigma/beta - 1;
@#else
r_star = mu_g^(1-gamma*(1-sigma))/beta - 1;
@#endif

model; //equation numbers refer to numbers in the working paper version
y=exp(z)*k(-1)^(1-alpha)*(g*l)^alpha; // Production technology (1)
z = rho_z*z(-1)+sigma_z*eps_z; // Transitory shock (2)
log(g) = (1-rho_g)*log(mu_g)+rho_g*log(g(-1))+sigma_g*eps_g; // Trend shock
@#if ghh == 1
u = (c-tau*l^nu)^(1-sigma)/(1-sigma); // GHH utility (3)
uc = (c - tau*l^nu)^(-sigma);
ul = -tau*nu*l^(nu-1)*(c - tau*l^nu)^(-sigma);
f = beta*g^(1-sigma);
@#else
u = (c^gamma*(1-l)^(1-gamma))^(1-sigma)/(1-sigma); // Cobb-Douglas utility (4)
uc = gamma*u/c*(1-sigma);
ul = -(1-gamma)*u/(1-l)*(1-sigma);
f = beta*g^(gamma*(1-sigma));
@#endif
c+g*k=y+(1-delta)*k(-1)-phi/2*(g*k/k(-1)-mu_g)^2*k(-1)-b(-1)+q*g*b; // Resource constraint (5)
1/q = 1+r_star+psi*(exp(b-b_star)-1); // Price of debt (6)
uc*(1+phi*(g*k/k(-1)-mu_g))*g=f*uc(+1)*(1-delta+(1-alpha)*y(+1)/k+phi/2*(g(+1)*k(+1)/k-mu_g)*(g(+1)*k(+1)/k+mu_g)); // FOC wrt to capital (10) with envelope condition plugged in
ul+uc*alpha*y/l=0; // Leisure-consumption arbitrage (11)
uc*g*q=f*uc(+1); // Euler equation (12)

//definition of auxilary variables to be plotted
tb_y = (b(-1)-g*q*b)/y; // Trade balance to GDP ratio, not logged as it can be negative
c_y = log(c/y); // Consumption to GDP ratio, logged to be in percent
i_y = log((g*k-(1-delta)*k(-1)+phi/2*(g*k/k(-1)-mu_g)^2*k(-1))/y); // Investment to GDP ratio, logged to be in percent
end;

initval;
q = 1/(1+r_star);
b = b_star;
z = 0;
g = mu_g;

c = 0.583095;
k = 4.02387;
y = 0.721195;
l = 0.321155;

@#if ghh == 1
u = (c-tau*l^nu)^(1-sigma)/(1-sigma);
uc = (c - tau*l^nu)^(-sigma);
ul = -tau*nu*l^(nu-1)*(c - tau*l)^(-sigma);
f = beta*g^(1-sigma);
@#else
u = (c^gamma*(1-l)^(1-gamma))^(1-sigma)/(1-sigma);
uc = gamma*u/c*(1-sigma);
ul = -(1-gamma)*u/(1-l)*(1-sigma);
f = beta*g^(gamma*(1-sigma));
@#endif

tb_y = (b-g*q*b)/y;
c_y = c/y;
i_y = (g*k-(1-delta)*k)/y;
end;

shocks;
var eps_g; stderr 1/sigma_g/100; // use a 1 percent shock
var eps_z; stderr 1/sigma_z/100; // use a 1 percent shock
end;

steady;

check;

// Plot impulse response functions (Figure 4)
stoch_simul(order=1) tb_y c_y i_y;
