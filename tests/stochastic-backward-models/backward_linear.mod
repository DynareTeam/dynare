var y gdp gdp_pot g pi P rs;
varexo e_y e_g e_p;

parameters alpha_1 alpha_2 beta_1 beta_2 theta gamma_1 gamma_2 pi_tar g_ss rr_bar;
alpha_1 = 0.9;
alpha_2 = -0.1;
beta_1 = 0.8;
beta_2 = 0.1;
theta = 0.9;
gamma_1 = 1.5;
gamma_2 = 0.5;
pi_tar = 2;
g_ss = 0.005;
rr_bar = 1;

model(linear);
gdp = gdp_pot + y;
y = alpha_1*y(-1) + alpha_2*(rs - pi - rr_bar) + e_y;
gdp_pot = gdp_pot(-1) + g;
g = (1 - theta)*g_ss + theta*g(-1) + e_g;
pi = (1 - beta_1)*pi_tar + beta_1*pi(-1) + beta_2*y + e_p;
P = pi/400 + P(-1);
rs = rr_bar + pi_tar + gamma_1*(pi(-1) - pi_tar) + gamma_2*y(-1);
end;

histval;
gdp_pot(0) = 1;
P(0) = 1;
g(0) = g_ss;
pi(0) = pi_tar;
end;

oo_.steadystate = NaN(7,1);

oo_ = simul_backward_model(M_.endo_histval, 10, options_, M_, oo_, zeros(11,3));

err1 = norm(abs(oo_.endo_simul([1 4 5 7],2:end) - repmat([0 0.005 2 3]',1,10)));
err2 = norm(abs(oo_.endo_simul([2 3 6],2:end) - repmat(linspace(1.005,1.05,10),3,1))); 
if err1 > 1e-14 || err2 > 1e-14;
   error('Error in backward_linear.mod');   
end;
