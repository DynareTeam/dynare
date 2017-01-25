// this program solves and simulates the model in
// "The Demand for Money during Hyperinflations under Rational Expectations: I" by T. Sargent, IER 1977
// this program mainly serves as the data generating process for the estimation of the model in sargent77ML.mod and sargent77Bayes.mod
// variables are defined as follows:
// x=p_t-p_{t-1}, p being the log of the price level
// mu=m_t-m_{t-1}, m being the log of money supply
// note that in contrast to the paper eta and epsilon have variance 1 (they are multiplied by the standard deviations)



var x mu a1 a2;
varexo epsilon eta;
parameters alpha lambda sig_eta sig_epsilon;
lambda=.5921;
alpha=-2.344;
sig_eta= .001;
sig_epsilon= .001;


// the model equations are taken from equation (27) on page 69 of the paper

model;
x=x(-1)-lambda*a1(-1)+(1/(lambda+alpha*(1-lambda)))*sig_epsilon*epsilon-(1/(lambda+alpha*(1-lambda)))*sig_eta*eta;
mu=(1-lambda)*x(-1)+lambda*mu(-1)-lambda*a2(-1)+(1+alpha*(1-lambda))/(lambda+alpha*(1-lambda))*sig_epsilon*epsilon-(1-lambda)/(lambda+alpha*(1-lambda))*sig_eta*eta;
a1=(1/(lambda+alpha*(1-lambda)))*sig_epsilon*epsilon-(1/(lambda+alpha*(1-lambda)))*sig_eta*eta;
a2=(1+alpha*(1-lambda))/(lambda+alpha*(1-lambda))*sig_epsilon*epsilon-(1-lambda)/(lambda+alpha*(1-lambda))*sig_eta*eta;
end;

steady;

shocks;

var eta;
stderr 1;
var epsilon;
stderr 1;
end;

stoch_simul(dr_algo=1,drop=0, order=1, periods=33, irf=0);

save data_hyperinfl.mat x mu;
