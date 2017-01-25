// Same model as sim_lead_lag_aux_vars.mod, but no auxiliary variables to deal with leads and lags

var c k z_forward z_backward;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1); // Resource constraint
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam); // Euler equation
z_backward=0.1*1+0.3*z_backward(-1)+0.3*z_backward(-2)+0.3*z_backward(-3)+(x(-4)-1);
z_forward=0.1*1+0.45*z_forward(+1)+0.45*z_forward(+2)+(x(+4)-1);
end;

initval;
c = 1.2;
k = 12;
x = 1; %set x steady state
end;

histval;
k(0) = 12;
x(0) = 1;
x(-1)=1.30; %set x(-1)
x(-2)=1.30; %set x(-2)
z_backward(0)=0.5; %set z_backward(0)
z_backward(-1)=0.4; %set z_backward(-1)
z_backward(-2)=0.9; %set z_backward(-2)
end;


shocks;
var x; %set x(2)
periods 2;
values 0.9;
end;

simul(periods=200,maxit=100);
base_results=load('sim_base_results.mat');
if max(abs(base_results.oo_.endo_simul(strmatch('c',base_results.M_.endo_names,'exact'),1+base_results.M_.maximum_endo_lag:end-base_results.M_.maximum_endo_lead) -...
    oo_.endo_simul(strmatch('c',M_.endo_names,'exact'),1+M_.maximum_lag:end-M_.maximum_lead)))>1e-8 || ...
    max(abs(base_results.oo_.endo_simul(strmatch('k',base_results.M_.endo_names,'exact'),1+base_results.M_.maximum_endo_lag:end-base_results.M_.maximum_endo_lead) -...
        oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),1+M_.maximum_lag:end-M_.maximum_lead)))>1e-8
    error('Autonomous system part is wrong')
end

clear base_results
base_results_aux_vars=load('sim_lead_lag_aux_vars_results.mat');

if max(abs(base_results_aux_vars.oo_.endo_simul(strmatch('x_lag_3',base_results_aux_vars.M_.endo_names,'exact'),1:end-base_results_aux_vars.M_.maximum_lead-1)' -...
    oo_.exo_simul(1:end-M_.maximum_lead-M_.maximum_lag,strmatch('x',M_.exo_names,'exact'))))>1e-8 
    error('Translation of aux vars is wrong')
end