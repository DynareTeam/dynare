periods 5000;

var c k mu_c b d in;
varexo e_d e_b; 

parameters R rho rho_b mu_b mu_d; 
R=1.05;
//rho=0.9;
rho = 0;
mu_b=30;
mu_d=5;
rho_b = 0;

model(linear);

 c+k = R*k(-1) + d;
 mu_c = b - c;
 mu_c=mu_c(+1);
 d= rho*d(-1)+ mu_d*(1-rho) + e_d;
 b=(1-rho_b)*mu_b+rho_b*b(-1)+e_b;
 in = k - k(-1);
 end;

//With a unit root, there exists no steady state.  Use the following trick.
//Supply ONE solution corresponding to the initial k that you named.

initval;
d=mu_d;
k=100;
c = (R-1)*k +d;
mu_c=mu_b-c;
b=mu_b;
end;

shocks;
var e_d;
stderr 1;
var e_b;
stderr 1;
end;

steady;
check;

stoch_simul(dr_algo=1, order=1, periods=500, irf=10);
save data_hall.mat c in;
