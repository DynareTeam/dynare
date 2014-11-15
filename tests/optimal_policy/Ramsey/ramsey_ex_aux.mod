/* Mod file tests the correctness of the Ramsey command when used together with a steady state file by
 * - checking whether the results coincide with the ones when used with an initval block
 * - checking whether between stoch_simul and ramsey_planner are consistent
 * The example is taken from Juillard, Michel (2011): User manual for optimal policy package, 
 * MONFISPOL FP7 project SSH-225149, Deliverable 1.1.2
*/

var pai, c, n, r, a, junk_endo_backward,junk_endo_forward,junk_expectation,junk_exo_backward,junk_exo_forward;
varexo u;
parameters beta, rho, epsilon, omega, phi, gamma;

beta=0.99;
gamma=3;
omega=17;
epsilon=8;
phi=1;
rho=0.95;

model;
a = rho*a(-1)+u;
1/c = beta*r/(c(+1)*pai(+1));
pai*(pai-1)/c = beta*pai(+1)*(pai(+1)-1)/c(+1)+epsilon*phi*n^(gamma+1)/omega -exp(a)*n*(epsilon-1)/(omega*c);
junk_endo_backward=c(+2);
junk_endo_forward=c(-2);
junk_expectation=EXPECTATION(-1)(c(+1));
junk_exo_backward=u(-2);
junk_exo_forward=u(+2);
exp(a)*n = c+(omega/2)*(pai-1)^2;
end;

initval;
r=1;
end;

steady_state_model;
a = 0;
pai = beta*r;
c = find_c(0.96,pai,beta,epsilon,phi,gamma,omega);
n = c+(omega/2)*(pai-1)^2;
junk_endo_backward=c;
junk_endo_forward=c;
junk_expectation=c;
junk_exo_backward=0;
junk_exo_forward=0;
end;

shocks;
var u; stderr 0.008;
end;

planner_objective(ln(c)-phi*((n^(1+gamma))/(1+gamma)));
ramsey_policy(planner_discount=0.99,order=1,instruments=(r));
