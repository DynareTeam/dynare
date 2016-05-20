// Example of optimal simple rule
var y inflation r dummy_var;
varexo y_ inf_;

parameters delta sigma alpha kappa gammax0 gammac0 gamma_y_ gamma_inf_;

delta =  0.44;
kappa =  0.18;
alpha =  0.48;
sigma = -0.06;


model(linear);
y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - inflation(+1)) + y_; 
inflation  =   alpha * inflation(-1) + (1-alpha) * inflation(+1) + kappa*y + inf_;
dummy_var=0.9*dummy_var(-1)+0.01*y;
r = gammax0*y(-1)+gammac0*inflation(-1)+gamma_y_*y_+gamma_inf_*inf_;
end;

shocks;
var y_;
stderr 0.63;
var inf_;
stderr 0.4;
end;

options_.nograph=1;
options_.nocorr=1;
options_.osr.tolf=1e-20;
osr_params gammax0 gammac0 gamma_y_ gamma_inf_;


optim_weights;
inflation 1;
y 1;
dummy_var 1;
end;


gammax0 = 0.2;
gammac0 = 1.5;
gamma_y_ = 8;
gamma_inf_ = 3;

osr;
%compute objective function manually
objective=oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('y',M_.endo_names,'exact'))+oo_.var(strmatch('inflation',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'))+oo_.var(strmatch('dummy_var',M_.endo_names,'exact'),strmatch('dummy_var',M_.endo_names,'exact'));
        
if abs(oo_.osr.objective_function-objective)>1e-8
    error('Objective Function is wrong')
end

%redo computation with covariance specified
optim_weights;
inflation 1;
y 1;
dummy_var 1;
y,inflation 0.5;
end;

osr;
%compute objective function manually
objective=oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('y',M_.endo_names,'exact'))+oo_.var(strmatch('inflation',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'))+oo_.var(strmatch('dummy_var',M_.endo_names,'exact'),strmatch('dummy_var',M_.endo_names,'exact'))+0.5*oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'));
if abs(oo_.osr.objective_function-objective)>1e-8
    error('Objective Function is wrong')
end

gammax0=1.35533;
gammac0=1.39664;
gamma_y_=16.6667;
gamma_inf_=9.13199;
        
%redo computation with double weight on one covariance 
optim_weights;
inflation 1;
y 1;
dummy_var 1;
y,inflation 1;
end;

osr;
%compute objective function manually
objective=oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('y',M_.endo_names,'exact'))+oo_.var(strmatch('inflation',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'))+oo_.var(strmatch('dummy_var',M_.endo_names,'exact'),strmatch('dummy_var',M_.endo_names,'exact'))+1*oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'));
if abs(oo_.osr.objective_function-objective)>1e-8
    error('Objective Function is wrong')
end
oo_covar_single=oo_;

%redo computation with single weight on both covariances

optim_weights;
inflation 1;
y 1;
dummy_var 1;
y,inflation 0.5;
inflation,y 0.5;
end;

osr;
%compute objective function manually
objective=oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('y',M_.endo_names,'exact'))+oo_.var(strmatch('inflation',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'))+oo_.var(strmatch('dummy_var',M_.endo_names,'exact'),strmatch('dummy_var',M_.endo_names,'exact'))+0.5*oo_.var(strmatch('y',M_.endo_names,'exact'),strmatch('inflation',M_.endo_names,'exact'))+0.5*oo_.var(strmatch('inflation',M_.endo_names,'exact'),strmatch('y',M_.endo_names,'exact'));
if abs(oo_.osr.objective_function-objective)>1e-8
    error('Objective Function is wrong')
end
if abs(oo_.osr.objective_function-oo_covar_single.osr.objective_function)>1e-8
    error('Objective Function is wrong')
end
if max(abs(cell2mat(struct2cell(oo_.osr.optim_params))-cell2mat(struct2cell(oo_covar_single.osr.optim_params))))>1e-5
    error('Parameters should be identical')
end
