function [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(xx,info,mh_conf_sig)
% Computes posterior mean, median, variance, HPD interval, deciles, and density from posterior draws.
%
% INPUTS
%    xx            [double]    Vector of posterior draws (or prior draws ;-)
%    info          [integer]   If equal to one the posterior density is estimated.
%    mh_config_sig [double]    Scalar between 0 and 1 specifying the size of the HPD interval.
%
%
% OUTPUTS
%    post_mean     [double]    Scalar, posterior mean.
%    post_median   [double]    Scalar, posterior median.
%    post_var      [double]    Scalar, posterior variance.
%    hpd_interval  [double]    Vector (1*2), Highest Probability Density interval 
%    post_deciles  [double]    Vector (9*1), deciles of the posterior distribution.
%    density       [double]    Matrix (n*2), non parametric estimate of the posterior density. First and second 
%                              columns are respectively abscissa and ordinate coordinates.  
%    
% SPECIAL REQUIREMENTS
%    Other matlab routines distributed with Dynare: mh_optimal_bandwidth.m 
%                                                   kernel_density_estimate.m.
%

% Copyright (C) 2005-2011 Dynare Team
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

xx = xx(:);
xx = sort(xx);


post_mean = mean(xx);
post_median = median(xx);
post_var = var(xx);

number_of_draws = length(xx);
hpd_draws = round((1-mh_conf_sig)*number_of_draws);

if hpd_draws>2
    kk = zeros(hpd_draws,1);
    jj = number_of_draws-hpd_draws;
    for ii = 1:hpd_draws
        kk(ii) = xx(jj)-xx(ii);
        jj = jj + 1;
    end
    [kmin,idx] = min(kk);
    hpd_interval = [xx(idx) xx(idx)+kmin];
else
    hpd_interval=NaN(1,2);    
end
if length(xx)>9
    post_deciles = xx([round(0.1*number_of_draws) ...
                       round(0.2*number_of_draws) ...
                       round(0.3*number_of_draws) ...
                       round(0.4*number_of_draws) ...
                       round(0.5*number_of_draws) ...
                       round(0.6*number_of_draws) ...
                       round(0.7*number_of_draws) ...
                       round(0.8*number_of_draws) ...
                       round(0.9*number_of_draws)]);
else
    post_deciles=NaN(9,1);
end
density = [];
if info
    number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
    if post_var > 1e-12
        bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
        kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.  
        optimal_bandwidth = mh_optimal_bandwidth(xx,number_of_draws,bandwidth,kernel_function);
        [density(:,1),density(:,2)] = kernel_density_estimate(xx,number_of_grid_points,...
                                                          number_of_draws,optimal_bandwidth,kernel_function);
    else
        density = NaN(number_of_grid_points,2);
    end
end