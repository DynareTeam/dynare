function  [raftery_lewis] = raftery_lewis(runs,q,r,s)
% function  raftery_lewis = raftery_lewis(runs,q,r,s)
% Computes the convergence diagnostics of Raftery and Lewis (1992), i.e. the 
% number of draws needed in MCMC to estimate the posterior cdf of the q-quantile 
% within an accuracy r with probability s
% 
% Inputs:
%   - draws [n_draws by n_var]      double matrix of draws from the sampler 
%   - q     [scalar]                quantile of the quantity of interest
%   - r     [scalar]                level of desired precision 
%   - s     [scalar]                probability associated with r
% 
% Output:
%   raftery_lewis   [structure]     containing the fields:
%   - M_burn    [n_draws by 1]      number of draws required for burn-in
%   - N_prec    [n_draws by 1]      number of draws required to achieve desired precision r
%   - k_thin    [n_draws by 1]      thinning required to get 1st order MC
%   - k_ind     [n_draws by 1]      thinning required to get independence
%   - I_stat    [n_draws by 1]      I-statistic of Raftery/Lewis (1992b)
%                                   measures increase in required
%                                   iterations due to dependence in chain
%   - N_min     [scalar]            # draws if the chain is white noise
%   - N_total   [n_draws by 1]      nburn + nprec
% 

% ---------------------------------------------------------------------
% NOTES:   Example values of q, r, s:
%     0.025, 0.005,  0.95 (for a long-tailed distribution)
%     0.025, 0.0125, 0.95 (for a short-tailed distribution);
%
%  - The result is quite sensitive to r, being proportional to the
%       inverse of r^2.
%  - For epsilon (closeness of probabilities to equilibrium values), 
%       Raftery/Lewis use 0.001 and argue that the results
%       are quite robust to changes in this value
%
% ---------------------------------------------------------------------
% REFERENCES: 
% Raftery, Adrien E./Lewis, Steven (1992a): "How many iterations in the Gibbs sampler?"
%   in: Bernardo/Berger/Dawid/Smith (eds.): Bayesian Statistics, Vol. 4, Clarendon Press: Oxford, 
%   pp. 763-773.   
% Raftery, Adrien E./Lewis, Steven (1992b): "Comment: One long run with diagnostics: 
%   Implementation strategies for Markov chain Monte Carlo." Statistical Science, 
%   7(4), pp. 493-497.  
%
% ----------------------------------------------------

% Copyright (C) 2016 Benjamin Born and Johannes Pfeifer
% Copyright (C) 2016 Dynare Team
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



[n_runs, n_vars] = size(runs);

raftery_lewis.M_burn=NaN(n_vars,1);
raftery_lewis.N_prec=NaN(n_vars,1);
raftery_lewis.k_thin=NaN(n_vars,1);
raftery_lewis.k_ind=NaN(n_vars,1);
raftery_lewis.I_stat=NaN(n_vars,1);
raftery_lewis.N_total=NaN(n_vars,1);


thinned_chain  = zeros(n_runs,1);
%quantities that can be precomputed as they are independent of variable
Phi   = norminv((s+1)/2); %note the missing ^{-1} at the Phi in equation top page 5, see RL (1995)
raftery_lewis.N_min  = fix(Phi^2*(1-q)*q/r^2+1);

for nv = 1:n_vars % big loop over variables    
    if q > 0 && q < 1
        work = (runs(:,nv) <= quantile(runs(:,nv),q));
    else
        error('Quantile must be between 0 and 1'); 
    end;
    
    k_thin_current_var = 1; 
    bic = 1; 
    epss = 0.001;
    % Find thinning factor for which first-order Markov Chain is preferred to second-order one
    while(bic > 0)
        thinned_chain=work(1:k_thin_current_var:n_runs,1); 
        [g2, bic] = first_vs_second_order_MC_test(thinned_chain);
        k_thin_current_var = k_thin_current_var+1;
    end;
    
    k_thin_current_var = k_thin_current_var-1; %undo last step
    
    %compute transition probabilities 
    transition_matrix = zeros(2,2);
    for i1 = 2:size(thinned_chain,1)
        transition_matrix(thinned_chain(i1-1)+1,thinned_chain(i1)+1) = transition_matrix(thinned_chain(i1-1)+1,thinned_chain(i1)+1)+1;
    end;
    alpha = transition_matrix(1,2)/(transition_matrix(1,1)+transition_matrix(1,2)); %prob of going from 1 to 2
    beta = transition_matrix(2,1)/(transition_matrix(2,1)+transition_matrix(2,2));  %prob of going from 2 to 1

    kmind=k_thin_current_var;
    [g2, bic]=independence_chain_test(thinned_chain);
    
    while(bic > 0)
        thinned_chain=work(1:kmind:n_runs,1); 
        [g2, bic] = independence_chain_test(thinned_chain);
        kmind = kmind+1;
    end;
    
    m_star  = log((alpha + beta)*epss/max(alpha,beta))/log(abs(1 - alpha - beta)); %equation bottom page 4
    raftery_lewis.M_burn(nv) = fix((m_star+1)*k_thin_current_var);
    n_star  = (2 - (alpha + beta))*alpha*beta*(Phi^2)/((alpha + beta)^3 * r^2); %equation top page 5
    raftery_lewis.N_prec(nv) = fix(n_star+1)*k_thin_current_var;
    raftery_lewis.I_stat(nv) = (raftery_lewis.M_burn(nv) + raftery_lewis.N_prec(nv))/raftery_lewis.N_min;
    raftery_lewis.k_ind(nv)  = max(fix(raftery_lewis.I_stat(nv)+1),kmind);
    raftery_lewis.k_thin(nv) = k_thin_current_var;
    raftery_lewis.N_total(nv)= raftery_lewis.M_burn(nv)+raftery_lewis.N_prec(nv);
end;

end

function [g2, bic] = first_vs_second_order_MC_test(d)
%conducts a test of first vs. second order Markov Chain via BIC criterion
n_obs=size(d,1);
g2 = 0; 
tran=zeros(2,2,2);
for t_iter=3:n_obs     % count state transitions
    tran(d(t_iter-2,1)+1,d(t_iter-1,1)+1,d(t_iter,1)+1)=tran(d(t_iter-2,1)+1,d(t_iter-1,1)+1,d(t_iter,1)+1)+1;
end;
% Compute the log likelihood ratio statistic for second-order MC vs first-order MC. G2 statistic of Bishop, Fienberg and Holland (1975)
for ind_1 = 1:2
    for ind_2 = 1:2
        for ind_3 = 1:2
            if tran(ind_1,ind_2,ind_3) ~= 0
                fitted = (tran(ind_1,ind_2,1) + tran(ind_1,ind_2,2))*(tran(1,ind_2,ind_3) + tran(2,ind_2,ind_3))/...
                    (tran(1,ind_2,1) + tran(1,ind_2,2) + tran(2,ind_2,1) + tran(2,ind_2,2));
                focus = tran(ind_1,ind_2,ind_3);
                g2 = g2 + log(focus/fitted)*focus;
            end
        end;       % end of for i3
    end;        % end of for i2
end;         % end of for i1
g2 = g2*2;
bic = g2 - log(n_obs-2)*2;

end


function [g2, bic] = independence_chain_test(d)
%conducts a test of independence Chain via BIC criterion
n_obs=size(d,1);
trans = zeros(2,2);
for ind_1 = 2:n_obs
    trans(d(ind_1-1)+1,d(ind_1)+1)=trans(d(ind_1-1)+1,d(ind_1)+1)+1;
end;
dcm1 = n_obs - 1;  
g2 = 0;
% Compute the log likelihood ratio statistic for second-order MC vs first-order MC. G2 statistic of Bishop, Fienberg and Holland (1975)
for ind_1 = 1:2
    for ind_2 = 1:2
        if trans(ind_1,ind_2) ~= 0
            fitted = ((trans(ind_1,1) + trans(ind_1,2))*(trans(1,ind_2) + trans(2,ind_2)))/dcm1; 
            focus = trans(ind_1,ind_2);
            g2 = g2 + log(focus/fitted)*focus;
        end;
    end;
end;
g2 = g2*2; 
bic = g2 - log(dcm1);
end
