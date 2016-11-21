function [results_vec, results_struct] = geweke_moments(draws,Dynareoptions)
%[results_vec, results_struct] = geweke_moments(draws,Dynareoptions)
% PURPOSE: computes Gewke's convergence diagnostics NSE and RNE 
%          (numerical std error and relative numerical efficiencies)

% INPUTS 
%   draws            [ndraws by 1 vector] 
%   Dynareoptions    [structure]
%  
% OUTPUTS 
%   results_vec
%   results_struct   [structure]  containing the following fields:
%          posteriormean= posterior parameter mean
%          posteriorstd = posterior standard deviation
%          nse_iid      = nse assuming no serial correlation for variable i
%          rne_iid      = rne assuming no serial correlation for variable i
%          nse_x        = nse using x% autocovariance tapered estimate
%          rne_x        = rne using x% autocovariance taper
% -----------------------------------------------------------------

%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2013 Dynare Team
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

% REFERENCES: Geweke (1992), `Evaluating the accuracy of sampling-based
% approaches to the calculation of posterior moments', in J.O. Berger,
% J.M. Bernardo, A.P. Dawid, and A.F.M. Smith (eds.) Proceedings of
% the Fourth Valencia International Meeting on Bayesian Statistics,
% pp. 169-194, Oxford University Press
% Geweke (1999): `Using simulation methods for Bayesian econometric models: 
% Inference, development and communication', Econometric Reviews, 18(1),
% 1-73
% -----------------------------------------------------------------

% written by: Johannes Pfeifer, 
% based on code by James P. LeSage, who in turn 
% drew on MATLAB programs written by Siddartha Chib 

  
ndraw = size(draws,1);
n_groups=100;
taper_steps=Dynareoptions.convergence.geweke.taper_steps;
results_vec=zeros(1,4+2*length(taper_steps));

ns = floor(ndraw/n_groups); %step_size
n_draws_used = ns*n_groups; %effective number of draws used after rounding down

window_means= zeros(n_groups,1);
window_uncentered_variances= zeros(n_groups,1);
for ig=1:n_groups;
    window_means(ig,1)=sum(draws((ig-1)*ns+1:ig*ns,1))/ns;
    window_uncentered_variances(ig,1)=sum(draws((ig-1)*ns+1:ig*ns,1).^2)/ns;        
end; %for ig
total_mean=mean(window_means);
total_variance=mean(window_uncentered_variances)-total_mean^2;

% save posterior means and std deviations to results_struct structure
results_vec(1,1)=total_mean;
results_vec(1,2)=sqrt(total_variance);
results_struct.posteriormean = total_mean;
results_struct.posteriorstd = results_vec(1,2);

% numerical standard error assuming no serial correlation
NSE=std(draws(1:n_draws_used,1),1)/sqrt(n_draws_used);
% save to results_struct structure
results_vec(1,3)=NSE;
results_vec(1,4)=total_variance/(n_draws_used*NSE^2);
results_struct.nse_iid = NSE;
results_struct.rne_iid = results_vec(1,4);

%get autocovariance of grouped means
centered_window_means=window_means-total_mean;
autocov_grouped_means=zeros(n_groups,1);
for lag=0:n_groups-1;
    autocov_grouped_means(lag+1)=centered_window_means(lag+1:n_groups,1)'*centered_window_means(1:n_groups-lag,1)/100;
end;

% numerical standard error with tapered autocovariance functions
for taper_index=1:length(taper_steps)
    taper=taper_steps(taper_index);
    taper_lags=(1:taper-1)';
    taper_lag_weight=1-taper_lags/taper;
    tapered_sum_of_covariances=autocov_grouped_means(1)+sum(2*taper_lag_weight.*autocov_grouped_means(1+taper_lags));
    NSE_taper=sqrt(tapered_sum_of_covariances/n_groups);
    % save results_struct in structure
    results_vec(1,4+taper_index*2-1)=NSE_taper;
    results_vec(1,4+taper_index*2)=total_variance/(n_draws_used*NSE_taper^2);

    eval(['results_struct.nse_taper_',num2str(taper),'= NSE_taper;']);
    eval(['results_struct.rne_taper_',num2str(taper),'= total_variance/(n_draws_used*NSE_taper^2);']);
end; % end of for mm loop

