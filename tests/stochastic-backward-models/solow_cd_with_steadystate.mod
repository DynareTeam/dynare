var Efficiency                  // $A$
    EfficiencyGrowth            // $X$
    Population                  // $L$
    PopulationGrowth            // $N$
    Output                      // $Y$
    PhysicalCapitalStock ;      // $K$

varexo e_x   // $\varepsilon_x$
       e_n ; // $\varepsilon_n$

parameters alpha                               // $\alpha$
	   delta                               // $\delta$
	   s                                   // $s$
           rho_x                               // $\rho_x$
           rho_n                               // $\rho_n$
           EfficiencyGrowth_ss                 // $X^{\star}$
           PopulationGrowth_ss ;               // $N^{\star}$

alpha = .33;
delta = .02;
s     = .20;
rho_x = .90;
rho_n = .95;
EfficiencyGrowth_ss = 1.00; // Do not change this calibration
PopulationGrowth_ss = 1.00; // Do not change this calibration

model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = PhysicalCapitalStock(-1)^alpha*(Efficiency*Population)^(1-alpha);
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

histval;
    Efficiency(0) = .5;
    EfficiencyGrowth(0) = 1.00;
    Population(0) = 1;
    PopulationGrowth(0) = 1.04;
    PhysicalCapitalStock(0) = 15;
end;

LongRunEfficiency = M_.endo_histval(1)*M_.endo_histval(2)^(rho_x/(1-rho_x));
LongRunPopulation = M_.endo_histval(3)*M_.endo_histval(4)^(rho_n/(1-rho_n));
LongRunEfficiencyGrowth = EfficiencyGrowth_ss; 
LongRunPopulationGrowth = PopulationGrowth_ss;
LongRunIntensiveCapitalStock = LongRunEfficiencyGrowth*LongRunPopulationGrowth*(s/(LongRunEfficiencyGrowth*LongRunPopulationGrowth-1+delta))^(1/(1-alpha));

precision = 1e-6;
T = 5*floor(log(precision)/log(max(rho_x, rho_n)));

oo_ = simul_backward_model(M_.endo_histval, T, options_, M_, oo_, zeros(T+1,2));

if abs(oo_.endo_simul(1,end)-LongRunEfficiency)>1e-10
    error('Wrong long run level!')
end

if abs(oo_.endo_simul(3,end)-LongRunPopulation)>1e-10
    error('Wrong long run level!')
end

IntensiveCapitalStock = oo_.endo_simul(6,1:end)./(oo_.endo_simul(1,1:end).*oo_.endo_simul(3,1:end));

if abs(IntensiveCapitalStock-LongRunIntensiveCapitalStock)>1e-6
    error('Wrong long run level!')
end
