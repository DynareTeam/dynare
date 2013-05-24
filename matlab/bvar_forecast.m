function bvar_forecast(nlags)
% function bvar_forecast(nlags)
% builds forecats for a bvar model
%
% INPUTS
%    nlags:     number of lags for the bvar
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2007-2012 Dynare Team
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

global options_ oo_ M_

if options_.forecast == 0
    error('bvar_forecast: you must specify "forecast" option')
end
[ny, nx, posterior, prior, forecast_data] = bvar_toolbox(nlags);

sims_no_shock = NaN(options_.forecast, ny, options_.bvar_replic);
sims_with_shocks = NaN(options_.forecast, ny, options_.bvar_replic);

S_inv_upper_chol = chol(inv(posterior.S));

% Option 'lower' of chol() not available in old versions of
% Matlab, so using transpose
XXi_lower_chol = chol(posterior.XXi)';

k = ny*nlags+nx;

% Declaration of the companion matrix
Companion_matrix = diag(ones(ny*(nlags-1),1),-ny);

% Number of explosive VAR models sampled
p = 0;
% Loop counter initialization
d = 0;
while d <= options_.bvar_replic
    
    Sigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);

    % Option 'lower' of chol() not available in old versions of
    % Matlab, so using transpose
    Sigma_lower_chol = chol(Sigma)';

    Phi = rand_matrix_normal(k, ny, posterior.PhiHat, Sigma_lower_chol, XXi_lower_chol);
    
    % All the eigenvalues of the companion matrix have to be on or inside the unit circle
    Companion_matrix(1:ny,:) = Phi(1:ny*nlags,:)'; 
    test = (abs(eig(Companion_matrix)));
    if any(test>1.0000000000001)
        p = p+1;
    end
    d = d+1;
    
    % Without shocks
    lags_data = forecast_data.initval;
    for t = 1:options_.forecast
        X = [ reshape(flipdim(lags_data, 1)', 1, ny*nlags) forecast_data.xdata(t, :) ];
        y = X * Phi;
        lags_data(1:end-1,:) = lags_data(2:end, :);
        lags_data(end,:) = y;
        sims_no_shock(t, :, d) = y;
    end
    
    % With shocks
    lags_data = forecast_data.initval;
    for t = 1:options_.forecast
        X = [ reshape(flipdim(lags_data, 1)', 1, ny*nlags) forecast_data.xdata(t, :) ];
        shock = (Sigma_lower_chol * randn(ny, 1))';
        y = X * Phi + shock;
        lags_data(1:end-1,:) = lags_data(2:end, :);
        lags_data(end,:) = y;
        sims_with_shocks(t, :, d) = y;
    end
end

if p > 0
    disp(' ')
    disp(['Some of the VAR models sampled from the posterior distribution'])
    disp(['were found to be explosive (' num2str(p/options_.bvar_replic) ' percent).'])
    disp(' ')
end

% Plot graphs
sims_no_shock_mean = mean(sims_no_shock, 3);

sort_idx = round((0.5 + [-options_.conf_sig, options_.conf_sig, 0]/2) * options_.bvar_replic);

sims_no_shock_sort = sort(sims_no_shock, 3);
sims_no_shock_down_conf = sims_no_shock_sort(:, :, sort_idx(1));
sims_no_shock_up_conf = sims_no_shock_sort(:, :, sort_idx(2));
sims_no_shock_median = sims_no_shock_sort(:, :, sort_idx(3));

sims_with_shocks_sort = sort(sims_with_shocks, 3);
sims_with_shocks_down_conf = sims_with_shocks_sort(:, :, sort_idx(1));
sims_with_shocks_up_conf = sims_with_shocks_sort(:, :, sort_idx(2));

OutputDirectoryName = CheckPath('graphs',M_.fname);

dyn_graph=dynare_graph_init(sprintf('BVAR forecasts (nlags = %d)', nlags), ny, {'b-' 'g-' 'g-' 'r-' 'r-'});

for i = 1:ny
    dyn_graph=dynare_graph(dyn_graph,[ sims_no_shock_median(:, i) ...
                   sims_no_shock_up_conf(:, i) sims_no_shock_down_conf(:, i) ...
                   sims_with_shocks_up_conf(:, i) sims_with_shocks_down_conf(:, i) ], ...
                 options_.varobs(i, :));
end

dyn_saveas(dyn_graph.fh,[OutputDirectoryName '/' M_.fname '_BVAR_forecast_',num2str(nlags)],options_)

% Compute RMSE

if ~isempty(forecast_data.realized_val)
    
    sq_err_cumul = zeros(1, ny);
    
    lags_data = forecast_data.initval;
    for t = 1:size(forecast_data.realized_val, 1)
        X = [ reshape(flipdim(lags_data, 1)', 1, ny*nlags) forecast_data.realized_xdata(t, :) ];
        y = X * posterior.PhiHat;
        lags_data(1:end-1,:) = lags_data(2:end, :);
        lags_data(end,:) = y;
        sq_err_cumul = sq_err_cumul + (y - forecast_data.realized_val(t, :)) .^ 2;
    end
    
    rmse = sqrt(sq_err_cumul / size(forecast_data.realized_val, 1));
    
    fprintf('RMSE of BVAR(%d):\n', nlags);
    
    for i = 1:size(options_.varobs, 1)
        fprintf('%s: %10.4f\n', options_.varobs(i, :), rmse(i));
    end 
end

% Store results

DirectoryName = [ M_.fname '/bvar_forecast' ];
if ~isdir(DirectoryName)
    if ~isdir(M_.fname)
        mkdir(M_.fname);
    end
    mkdir(DirectoryName);
end
save([ DirectoryName '/simulations.mat'], 'sims_no_shock', 'sims_with_shocks');

for i = 1:size(options_.varobs, 1)
    name = options_.varobs(i, :);

    sims = squeeze(sims_with_shocks(:,i,:));
    eval(['oo_.bvar.forecast.with_shocks.Mean.' name ' = mean(sims, 2);']);
    eval(['oo_.bvar.forecast.with_shocks.Median.' name ' = median(sims, 2);']);
    eval(['oo_.bvar.forecast.with_shocks.Var.' name ' = var(sims, 0, 2);']);
    eval(['oo_.bvar.forecast.with_shocks.HPDsup.' name ' = sims_with_shocks_up_conf(:,i);']);
    eval(['oo_.bvar.forecast.with_shocks.HPDinf.' name ' = sims_with_shocks_down_conf(:,i);']);

    sims = squeeze(sims_no_shock(:,i,:));
    eval(['oo_.bvar.forecast.no_shock.Mean.' name ' = sims_no_shock_mean(:,i);']);
    eval(['oo_.bvar.forecast.no_shock.Median.' name ' = sims_no_shock_median(:,i);']);
    eval(['oo_.bvar.forecast.no_shock.Var.' name ' = var(sims, 0, 2);']);
    eval(['oo_.bvar.forecast.no_shock.HPDsup.' name ' = sims_no_shock_up_conf(:,i);']);
    eval(['oo_.bvar.forecast.no_shock.HPDinf.' name ' = sims_no_shock_down_conf(:,i);']);

    if exist('rmse')
        eval(['oo_.bvar.forecast.rmse.' name ' = rmse(i);']);
    end
end
