/*
 * Example 1 from F. Collard (2001): "Stochastic simulations with DYNARE:
 * A practical guide" (see "guide.pdf" in the documentation directory).
 */

/*
 * Copyright (C) 2001-2010 Dynare Team
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


var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((b*c)/(b(+1)*c(+1)))
    *(b(+1)*alpha*y(+1)+(1-delta)*k));
y = a*(k(-1)^alpha)*(h^(1-alpha));
k = b*(y-c)+(1-delta)*k(-1);
log(a) = rho*log(a(-1))+tau*log(b(-1)) + e;
log(b) = tau*log(a(-1))+rho*log(b(-1)) + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 1;
b = 1;
end;
resid(1);
shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul(loglinear,order=1);

load results_exp;
if max(max(abs(oo_.dr.ghx-oo_exp.dr.ghx)))>1e-10
    error('Option loglinear wrong, ghx not equal')
end
if max(max(abs(oo_.dr.ghu-oo_exp.dr.ghu)))>1e-10
    error('Option loglinear wrong, ghu not equal')
end
if max(max(abs(oo_.irfs.y_e-oo_exp.irfs.y_e)))>1e-10
    error('Option loglinear wrong, IRFs not equal')
end
if max(max(abs(oo_.irfs.y_u-oo_exp.irfs.y_u)))>1e-10
    error('Option loglinear wrong, ghu not equal')
end
if max(max(abs(oo_.mean-oo_exp.mean)))>1e-10
    error('Option loglinear wrong, mean not equal')
end
if max(max(abs(oo_.dr.ys-oo_exp.dr.ys)))>1e-10
    error('Option loglinear wrong, ys not equal')
end
if max(max(abs(oo_.steady_state-oo_exp.steady_state)))>1e-10
    error('Option loglinear wrong, steady_state not equal')
end

for ii=1:length(oo_.gamma_y)
    if max(max(abs(oo_.gamma_y{ii,1}-oo_exp.gamma_y{ii,1})))>1e-10
    error('Option loglinear wrong, moments not equal')
    end
end

for ii=1:length(oo_.autocorr)
    if max(max(abs(oo_.autocorr{1,ii}-oo_exp.autocorr{1,ii})))>1e-10
    error('Option loglinear wrong, moments not equal')
    end
end
stoch_simul(loglinear,order=1,periods=100000);
if abs(mean(y)-0.0776)>0.02
    error('Simulations are wrong')
end