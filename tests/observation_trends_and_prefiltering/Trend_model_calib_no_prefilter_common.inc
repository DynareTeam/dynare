var Y_obs P_obs junk1 junk2;
varexo e_y e_p eps_junk;

parameters rho_y rho_p  g_y g_p const_y const_p sigma_y sigma_p;

rho_y=0.5;
rho_p=0.5;
g_y=0.0001;
g_p=-0.0001;
const_y=2;
const_p=2;
sigma_y=0.001;
sigma_p=0.001;


model;
Y_obs = (1-rho_y)*const_y + rho_y*Y_obs(-1)+sigma_y*e_y;
P_obs = (1-rho_p)*const_p + rho_p*P_obs(-1)+sigma_p*e_p;
junk1 = 0.9*junk1(+1);
junk2 = 0.9*junk2(-1)+eps_junk;
end;

steady_state_model;
Y_obs = const_y;
P_obs = const_p;
junk1=0;
junk2=0;
end;

shocks;
var e_p; stderr 1;
var e_y; stderr 1;
var eps_junk; stderr 1;
end;

steady;

varobs P_obs Y_obs junk2;

observation_trends;
P_obs (g_p);
Y_obs (g_y);
end;
