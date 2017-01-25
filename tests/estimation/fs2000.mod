// See fs2000.mod in the examples/ directory for details on the model

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
end;

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
end;


shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;

estimated_params;
alp, beta_pdf, 0.356, 0.02;
bet, beta_pdf, 0.993, 0.002;
gam, normal_pdf, 0.0085, 0.003;
mst, normal_pdf, 1.0002, 0.007;
rho, beta_pdf, 0.129, 0.223;
psi, beta_pdf, 0.65, 0.05;
del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
end;

varobs gp_obs gy_obs;

options_.solve_tolf = 1e-12;

estimation(order=1,datafile=fsdat_simul,nobs=192,loglinear,mh_replic=3000,mh_nblocks=1,mh_jscale=0.8,moments_varendo,selected_variables_only,contemporaneous_correlation,smoother,forecast=8,
        geweke_interval = [0.19 0.49],
        taper_steps = [4 7 15],
        raftery_lewis_diagnostics,
        raftery_lewis_qrs=[0.025 0.01 0.95],
        bayesian_irf,posterior_nograph
        ) y m;

if ~isequal(options_.convergence.geweke.taper_steps,[4 7 15]') || ~isequal(options_.convergence.geweke.geweke_interval,[0.19 0.49])
    error('Interface for Geweke diagnostics not working')
end
        
if ~isequal(options_.convergence.rafterylewis.qrs,[0.025 0.01 0.95]) || ~isequal(options_.convergence.rafterylewis.indicator,1)
    error('Interface for Raftery/Lewis diagnostics not working')
end

%test load_mh_file option
options_.bayesian_irf=0;
options_.smoother=0;
options_.moments_varendo=0;
options_.forecast=0;   
copyfile([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat'],[M_.dname '_mh1_blck1.mat'])
estimation(mode_compute=0,mode_file=fs2000_mode,order=1, datafile=fsdat_simul, nobs=192, loglinear, mh_replic=1500, mh_nblocks=1, mh_jscale=0.8);
hh=eye(size(bayestopt_.name,1));
save('fs2000_mode.mat','hh','-append')
Laplace = oo_.MarginalDensity.LaplaceApproximation; %save Laplace approximation which will be overwritten with set hh otherwise
estimation(mode_compute=0,mode_file=fs2000_mode,order=1, datafile=fsdat_simul, nobs=192, loglinear, mh_replic=1500, mh_nblocks=1, mh_jscale=10,load_mh_file);

temp1=load([M_.dname '_mh1_blck1.mat']);
temp2=load([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat']);

if ~isoctave
    if max(max(abs(temp1.x2-temp2.x2)))>1e-10
        error('Adding draws did not result in the same chain')
    end
end
        
save('fs2000_result.mat','oo_')        
options_.load_results_after_load_mh=1;
estimation(mode_compute=0,mode_file=fs2000_mode,order=1, datafile=fsdat_simul, nobs=192, loglinear, mh_replic=0, mh_nblocks=1, mh_jscale=10,load_mh_file,smoother) gy_obs gp_obs;
oo_.MarginalDensity.LaplaceApproximation = Laplace; %reset correct Laplace
