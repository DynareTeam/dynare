var Efficiency                  // $A$
    EfficiencyGrowth            // $X$
    Population                  // $L$
    PopulationGrowth            // $N$
    Output                      // $Y$
    PhysicalCapitalStock ;      // $K$

varexo e_x   // $\varepsilon_x$
       e_n ; // $\varepsilon_n$

parameters alpha                               // $\alpha$
           epsilon                             // $\varepsilon$
	   delta                               // $\delta$
	   s                                   // $s$
           rho_x                               // $\rho_x$
           rho_n                               // $\rho_n$
           EfficiencyGrowth_ss                 // $X^{\star}$
           PopulationGrowth_ss ;               // $N^{\star}$

alpha   = .33;
epsilon = .70;
delta   = .02;
s       = .20;
rho_x   = .90;
rho_n   = .95;
EfficiencyGrowth_ss = 1.02;
PopulationGrowth_ss = 1.02;

model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = (alpha*PhysicalCapitalStock(-1)^((epsilon-1)/epsilon)+(1-alpha)*(Efficiency*Population)^((epsilon-1)/epsilon))^(epsilon/(epsilon-1));
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

histval;
    Efficiency(0) = 1;
    EfficiencyGrowth(0) = 1.02;
    Population(0) = 1;
    PopulationGrowth(0) = 1.02;
    PhysicalCapitalStock(0) = 1;
end;

shocks;
    var e_x = 0.005;
    var e_n = 0.001;
end;

oo_ = simul_backward_nonlinear_model([], 5000, options_, M_, oo_);
