// example 1 from Collard's guide to Dynare
var y, c, k, a, h, b;
varexo e,u;

parameters beta, rho, alpha, delta, theta, psi, tau, phi;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 5;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
//var e, u = phi*0.009*0.009;
end;

varobs a y;


@#for order in [1,2]
options_.order=@{order};
stoch_simul(conditional_variance_decomposition = 100,irf=0) a y k;
if max(abs(sum(oo_.variance_decomposition,2)-100))>1e-6
    error(['Variance decomposition at order ',num2str(options_.order),' does not work'])
end

stoch_simul(conditional_variance_decomposition = [1 2 3 5 10 100],irf=0) a y k;
if max(max(abs(sum(oo_.conditional_variance_decomposition,3)-1)))>1e-6
    error(['Conditional variance decomposition at order ',num2str(options_.order),' does not work'])
end

shocks;
var y; stderr 0.01;
var a; stderr 0.009;
end;

stoch_simul(conditional_variance_decomposition = [1 2 3 5 10 200],irf=0) a y k;
if max(max(abs(sum(oo_.conditional_variance_decomposition,3)-1)))>1e-6
    error(['Conditional variance decomposition at order ',num2str(options_.order),' does not work'])
end
if max(max(abs(sum(oo_.conditional_variance_decomposition_ME,3)-1)))>1e-6
    error(['Conditional variance decomposition at order ',num2str(options_.order),' does not work'])
end

nvar = size(var_list_,1);
SubsetOfVariables=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list_(i,:),M_.endo_names,'exact');
    SubsetOfVariables(i) = i_tmp;
end

[observable_pos,index_observables,index_subset]=intersect(SubsetOfVariables,options_.varobs_id,'stable');
y_pos=strmatch('y',var_list_,'exact');
y_pos_varobs=strmatch('y',options_.varobs,'exact');
a_pos_varobs=strmatch('a',options_.varobs,'exact');

if (oo_.conditional_variance_decomposition_ME(index_observables(y_pos_varobs),end,end)+oo_.var(y_pos,y_pos)/(oo_.var(y_pos,y_pos)+M_.H(y_pos_varobs,y_pos_varobs)))-1>1e-5...
        || abs(oo_.conditional_variance_decomposition_ME(index_observables(a_pos_varobs),1,end)-0.5)>1e-5    
    error(['Conditional variance decomposition at order ',num2str(options_.order),' with ME does not work'])
end

if (oo_.variance_decomposition_ME(index_observables(y_pos_varobs),end)/100+oo_.var(y_pos,y_pos)/(oo_.var(y_pos,y_pos)+M_.H(y_pos_varobs,y_pos_varobs)))-1>1e-5   
    error(['Unconditional variance decomposition at order ',num2str(options_.order),' with ME does not work'])
end

shocks;
var y; stderr 0;
end;

@#endfor

%% do simulated moments
@#for order in [1,2]
options_.order=@{order};

shocks;
var y; stderr 0.01;
var a; stderr 0.009;
end;

stoch_simul(irf=0,periods=1000000) a y k;
if max(abs(sum(oo_.variance_decomposition,2)-100))>2
    error(['Variance decomposition at order ',num2str(options_.order),' does not work'])
end

[observable_pos,index_observables,index_subset]=intersect(SubsetOfVariables,options_.varobs_id,'stable');
y_pos=strmatch('y',var_list_,'exact');
y_pos_varobs=strmatch('y',options_.varobs,'exact');
a_pos_varobs=strmatch('a',options_.varobs,'exact');

if (oo_.variance_decomposition_ME(index_observables(y_pos_varobs),end)/100+oo_.var(y_pos,y_pos)/(oo_.var(y_pos,y_pos)+M_.H(y_pos_varobs,y_pos_varobs)))-1>1e-5   
    error(['Unconditional variance decomposition at order ',num2str(options_.order),' with ME does not work'])
end

shocks;
var y; stderr 0;
end;

@#endfor