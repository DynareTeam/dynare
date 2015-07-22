var y c k m n R pi z u;
varexo	e sigma;	 
// sigma stands for phi in the eq 2.37 p.69

parameters alphha betta delta gamm phi1 eta a b rho  phi2 Psi thetass;  
//phi1 stands for capital phi in eq.2.68 and 2.69
//phi2 stands for lowercase phi in eq. 2.66

alphha = 0.36;
betta = 0.989; 
gamm = 0.5;
delta = 0.019;
phi1 = 2;
phi2 = 0;
eta = 1;
a = 0.95;
b = 2.56;
rho = 0.95;
Psi = 1.47630583;
thetass = 1.0125;

model;

(a*exp(c)^(1-b)+(1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))*a*exp(c)^(-b) = (a*exp(c)^(1-b)+(1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))*(1-a)*exp(m)^(-b)+betta*(a*exp(c(+1))^(1-b)+(1-a)*exp(m(+1))^(1-b))^((b-phi1)/(1-b))*a*exp(c(+1))^(-b)/(1+pi(+1));

Psi*(1-exp(n))^(-eta)/(a*exp(c)^(-b)*(a*exp(c)^(1-b) + (1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))) = (1-alphha)*exp(y)/exp(n);

(a*exp(c)^(1-b)+(1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))*a*exp(c)^(-b) = betta*exp(R(+1))*(a*exp(c(+1))^(1-b)+(1-a)*exp(m(+1))^(1-b))^((b-phi1)/(1-b))*a*exp(c(+1))^(-b);

exp(R) = alphha*exp(y)/exp(k(-1)) + 1-delta;

exp(k) = (1-delta)*exp(k(-1))+exp(y)-exp(c);

exp(y) = exp(z)*exp(k(-1))^alphha*exp(n)^(1-alphha);

exp(m) = exp(m(-1))*(u+thetass)/(1+pi);

z = rho*z(-1) + e;

u = gamm*u(-1) + phi2*z(-1) + sigma;

end;

shocks;
var e; stderr 0.007;
var sigma;stderr 0.0089;
end;


steady;

