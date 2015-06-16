% Real Business Cycle Model

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k l x z g ysim xsim;

varexo eps_z eps_g;
parameters gn gz betta delta psi sigma theta rho_z eps_z_sigma rho_g eps_g_sigma q12 zbar gbar;


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

gn    = 0.00243918275778010;
gz    = 0.00499789993972673;
betta = 0.9722^(1/4)/(1+gz); 
delta = 1-(1-0.0464)^(1/4);
psi   = 2.24;
sigma = 1.000001;
theta = 0.35;

zbar  = 0.0023;  
gbar  = -0.0382;

rho_z = 0.8;
sig_z = 0.0086;
rho_g = 0.8;
sig_g = 0.0248;

q12 = -0.0002;



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------


model;
% Intratemporal Optimality Condition
(psi*c)/(1-l) = (1-theta)*(k(-1)^theta)*(l^(-theta))*(exp(z)^(1-theta)); 
% Intertemporal Optimality Condition
(c^(-sigma))*((1-l)^(psi*(1-sigma))) = betta*((c(+1))^(-sigma))*((1-l(+1))^(psi*(1-sigma)))*(theta*(k^(theta-1))*((exp(z(+1))*l(+1))^(1-theta))+(1-delta));
% Aggregate Resource Constraint
y = c + exp(g) + x;
% Capital Accumulation Law
(1+gz)*(1+gn)*k = (1-delta)*k(-1) + x;
% Production Function
y = (k(-1)^theta)*((exp(z)*l)^(1-theta));

% Stochastic Processes
z = zbar + rho_z*z(-1) + eps_z_sigma*eps_z;
g = gbar + rho_g*g(-1) + eps_g_sigma*eps_g;

ysim = y;
xsim = x;

end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------


shocks;
var eps_z ; stderr 1;
var eps_g ; stderr 1;
end;


initval;
z = 0;
k = (betta*theta)^(1/(1-theta))*exp(zbar/(1-rho_z));
c = (1-betta*theta)/(betta*theta)*k;
y = exp(zbar/(1-rho_z))^(1-theta)*k^theta;
l = k/(((1-betta*(1-delta))/(betta*theta*(exp(zbar/(1-rho_z))^(1-theta))))^(1/(theta-1)));
end;


estimated_params;
eps_z_sigma,    0.0001, 0,10;
eps_g_sigma,    0.0001, 0,10;
% betta, 0.9722^(1/4)/(1+gz),0,1;
%corr eps_z, eps_g,   0.0001, -1,1;
end;

varobs ysim xsim;

identification(advanced=1);
