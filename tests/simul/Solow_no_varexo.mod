/*
 * This file studies the transition to steady state in a simple Solow model with Cobb-Douglass production function
 *
 * This implementation was written by Johannes Pfeifer. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2014 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */



var c k y invest;
predetermined_variables k;
parameters s alpha delta n g;

s=0.2;
alpha=0.3;
delta=0.1;
n=0.01;
g=0.02;
g_initial=0.02;

model;
(1+n+g+n*g)*k(+1)=(1-delta)*k+invest;
invest+c=y;
c=(1-s)*y;
y=k^alpha;
end;

%set initial condition to 90% of steady state value
initval;
k=0.9*((delta+n+g+n*g)/s)^(1/(alpha-1));
y=k^alpha;
c=(1-s)*y;
invest=y-c;
end;

%set terminal condition to steady state value
endval;
k=((delta+n+g+n*g)/s)^(1/(alpha-1));
y=k^alpha;
c=(1-s)*y;
invest=y-c;
end;
resid(1);
steady;

simul(periods=200);
rplot k;
