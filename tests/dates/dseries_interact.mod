%% Mod-file tests interaction between estimation and shock_decomposition when dseries is used or not
var hh nn log_nn;
varexo eps_a;
parameters alfa nbar lambda betta rho_a std_a nn_init;

nn_init = -0.1;
alfa = 0.05;
lambda = 0.054;
betta = 0.99;
nbar = 1;
rho_a = 0;
std_a = 1;


model(linear);

hh = - alfa * nn + betta * ( hh(+1) + 0  *  eps_a(+1) ) + eps_a;

log_nn = log_nn(-1) + hh * lambda / (1-lambda);

log_nn = ln(nbar) + nn;

end;

steady_state_model;
log_nn = log(nbar);
nn = 0;
hh = 0;
end;

shocks;
var eps_a; stderr 1;
end;

estimated_params;
alfa, beta_pdf,   0.1, 0.05; 
std_a,   inv_gamma_pdf, 0.05, 1;
end;

varobs log_nn;

if ~isoctave() && ~matlab_ver_less_than('8.4')
   websave('data_uav.xlsx','http://www.dynare.org/Datasets/data_uav.xlsx', weboptions('Timeout', 30))
else
   urlwrite('http://www.dynare.org/Datasets/data_uav.xlsx','data_uav.xlsx')
end

%reading Excel sheet from column A on creates quarterly dseries starting in
%1950
estimation(first_obs=2,datafile=data_uav, xls_sheet=Tabelle1, xls_range=a1:b54, mh_replic=2, mh_nblocks=1, mh_jscale=1.1, mh_drop=0.8, plot_priors=0, smoother) log_nn nn hh ;
shock_decomposition( parameter_set=posterior_median ) nn hh;

%reading Excel sheet from column B on creates annual dseries starting with 1
estimation(first_obs=2,datafile=data_uav, xls_sheet=Tabelle1, xls_range=b1:b54, mh_replic=2, mh_nblocks=1, mh_jscale=1.1, mh_drop=0.8, plot_priors=0, smoother) log_nn nn hh ;
shock_decomposition( parameter_set=posterior_median ) nn hh;

delete('data_uav.xlsx')
