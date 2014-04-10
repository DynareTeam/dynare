// Example that triggers homotopy in perfect foresight simulation.

var Consumption, Capital, LoggedProductivity;

varexo LoggedProductivityInnovation;

parameters beta, alpha, delta, rho;

beta = .985;
alpha = 1/3;
delta = alpha/10;
rho = .9;

model;

[name='Euler equation'] // This is an equation tag!
1/Consumption = beta/Consumption(1)*(alpha*exp(LoggedProductivity(1))*Capital^(alpha-1)+1-delta);

[name='Physical capital stock law of motion']
Capital = exp(LoggedProductivity)*Capital(-1)^alpha+(1-delta)*Capital(-1)-Consumption;

[name='Logged productivity law of motion']
LoggedProductivity = rho*LoggedProductivity(-1)+LoggedProductivityInnovation;

end;

steady_state_model;
  LoggedProductivity = LoggedProductivityInnovation/(1-rho);
  Capital = (exp(LoggedProductivity)*alpha/(1/beta-1+delta))^(1/(1-alpha));
  Consumption = exp(LoggedProductivity)*Capital^alpha-delta*Capital;
end;

set_time(1Q1);

initval;
  LoggedProductivityInnovation = 0;
  LoggedProductivity = 10;
  Capital = 1;
end;

endval;
  LoggedProductivityInnovation = 0;
end;

steady;

simul(periods=200);

plot(Simulated_time_series.Capital(1Q1:25Q4));
