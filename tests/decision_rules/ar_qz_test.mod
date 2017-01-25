var white_noise ar1 junk;
varexo e;

parameters phi;

phi=1;

model;
white_noise=e;
ar1=phi*ar1(-1)+e;
junk=0.9*junk(+1);
end;

shocks;
var e = 1;
end;

options_.qz_criterium=1+1e-6;
stoch_simul(order=1);

options_.qz_criterium=1-1e-6;
error_indicator=0;
try 
    info=stoch_simul(var_list_)
    error_indicator=1
catch
    
end
if error_indicator
    error('qz_criterion did not work')
end