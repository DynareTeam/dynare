var y_backward dummy_var;

parameters rho_1 rho_2;

rho_1=0.2;
rho_2=0.1;

// Equilibrium conditions
model;
y_backward=rho_1*y_backward(-1)+rho_2*y_backward(-2);
dummy_var=0.9*dummy_var(+1);
end;

// Set starting value for solver
initval;
y_backward=1;
end;
steady;

// Set initial conditions for state variables

histval;
y_backward(0)=1;
y_backward(-1)=2;
end;


// Check the Blanchard-Kahn conditions
check;

// Deterministic simulation of the model for 200 periods
options_.solve_tolf=1e-12;
simul(periods=100);

// Display the path of consumption and capital
rplot y_backward;

junk=zeros(1,options_.periods+M_.maximum_lag);
junk(1)=2;
junk(2)=1;

for ii=3:options_.periods+M_.maximum_lag+M_.maximum_lead
    junk(ii)=M_.params(strmatch('rho_1',M_.param_names,'exact'))*junk(ii-1)+M_.params(strmatch('rho_2',M_.param_names,'exact'))*junk(ii-2);
end

if max(abs(junk(M_.maximum_lag+1:end)-oo_.endo_simul(strmatch('y_backward',M_.endo_names,'exact'),1:end-M_.maximum_lead)))>1e-10
    error('Solution of purely backwards model not correct')
end
        
