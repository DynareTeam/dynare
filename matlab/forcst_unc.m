function forcst_unc(y0,var_list)
% function [mean,intval1,intval2]=forcst_unc(y0,var_list)
% computes forecasts with parameter uncertainty
%
% INPUTS
%   y0: matrix of initial values
%   var_list: list of variables to be forecasted
%
% OUTPUTS
%   none
%
% ALGORITHM
%   uses antithetic draws for the shocks
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2013 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_ oo_ estim_params_ bayestopt_

% setting up estim_params_
[xparam1,estim_params_,bayestopt_,lb,ub] = set_prior(estim_params_,M_);

options_.TeX = 0;
options_.nograph = 0;
plot_priors(bayestopt_,M_,options_);

% workspace initialization
if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr,:);
end
n = size(var_list,1);

periods = options_.forecast;
exo_nbr = M_.exo_nbr;
replic = options_.replic;
order = options_.order;
maximum_lag = M_.maximum_lag;
%  params = prior_draw(1);
params = rndprior(bayestopt_);
set_parameters(params);
% eliminate shocks with 0 variance
i_exo_var = setdiff([1:exo_nbr],find(diag(M_.Sigma_e) == 0));
nx = length(i_exo_var);

ex0 = zeros(periods,exo_nbr);
yf1 = zeros(periods+M_.maximum_lag,n,replic);

% loops on parameter values
m1 = 0;
m2 = 0;
for i=1:replic
    % draw parameter values from the prior
    % params = prior_draw(0);
    params = rndprior(bayestopt_);
    set_parameters(params);
    % solve the model
    [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
    % discard problematic cases
    if info
        continue
    end
    % compute forecast with zero shocks
    m1 = m1+1;
    yf1(:,:,m1) = simult_(y0,dr,ex0,order)';
    % compute forecast with antithetic shocks
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
    ex(:,i_exo_var) = randn(periods,nx)*chol_S;
    m2 = m2+1;
    yf2(:,:,m2) = simult_(y0,dr,ex,order)';
    m2 = m2+1;
    yf2(:,:,m2) = simult_(y0,dr,-ex,order)';
end

oo_.forecast.accept_rate = (replic-m1)/replic;

if options_.noprint == 0 && m1 < replic
    disp(' ')
    disp(' ')
    disp('FORECASTING WITH PARAMETER UNCERTAINTY:')
    disp(sprintf(['The model  couldn''t be solved for %f%% of the parameter' ...
                  ' values'],100*oo_.forecast.accept_rate))
    disp(' ')
    disp(' ')
end

% compute moments
yf1 = yf1(:,:,1:m1);
yf2 = yf2(:,:,1:m2);

yf_mean = mean(yf1,3);

yf1 = sort(yf1,3);
yf2 = sort(yf2,3);

sig_lev = options_.conf_sig;
k1 = round((0.5+[-sig_lev, sig_lev]/2)*replic);
% make sure that lower bound is at least the first element
if k1(2) == 0
    k1(2) = 1;
end
k2 = round((1+[-sig_lev, sig_lev])*replic);
% make sure that lower bound is at least the first element
if k2(2) == 0
    k2(2) = 1;
end

% compute shock uncertainty around forecast with mean prior
set_parameters(bayestopt_.p1);
[dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
[yf3,yf3_intv] = forcst(dr,y0,periods,var_list);
yf3_1 = yf3'-[zeros(maximum_lag,n); yf3_intv];
yf3_2 = yf3'+[zeros(maximum_lag,n); yf3_intv];

% graphs
OutputDirectoryName = CheckPath('graphs',M_.fname);

dyn_graph=dynare_graph_init('Forecasts type I',n,{'b-' 'g-' 'g-' 'r-' 'r-'});
for i=1:n
    dynare_graph(dyn_graph,[yf_mean(:,i) squeeze(yf1(:,i,k1)) squeeze(yf2(:,i,k2))],...
                 var_list(i,:));
end
dyn_saveas(dyn_graph.fh,[OutputDirectoryName '/' M_.fname '_forecast_param_uncert_',num2str(nlags)],options_)

dyn_graph=dynare_graph_init('Forecasts type II',n,{'b-' 'k-' 'k-' 'r-' 'r-'});
for i=1:n
    dynare_graph(dyn_graph,[yf_mean(:,i) yf3_1(:,i) yf3_2(:,i) squeeze(yf2(:,i,k2))],...
                 var_list(i,:));
end
dyn_saveas(dyn_graph.fh,[OutputDirectoryName '/' M_.fname '_forecast_param_shock_uncert_',num2str(nlags)],options_)


% saving results
save_results(yf_mean,'oo_.forecast.Mean.',var_list);
save_results(yf1(:,:,k1(1)),'oo_.forecast.HPDinf.',var_list);
save_results(yf1(:,:,k1(2)),'oo_.forecast.HPDsup.',var_list);
save_results(yf2(:,:,k2(1)),'oo_.forecast.HPDTotalinf.',var_list);
save_results(yf2(:,:,k2(2)),'oo_.forecast.HPDTotalsup.',var_list);