function [options_, oo_]=ms_forecast(M_, options_, oo_)
% function [options_, oo_]=ms_forecast(M_, options_, oo_)
% Markov-switching SBVAR: Forecast
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2016 Dynare Team
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

disp('MS-SBVAR Forecasts');
options_ = set_file_tags(options_);
[options_, oo_] = set_ms_estimation_file(options_.ms.file_tag, options_, oo_);
clean_ms_forecast_files(options_.ms.output_file_tag);
forecastdir = [options_.ms.output_file_tag filesep 'Forecast'];
create_dir(forecastdir);

% setup command line options
opt = ['-forecast -nodate -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -do ' forecastdir];
opt = [opt ' -ft ' options_.ms.file_tag];
opt = [opt ' -fto ' options_.ms.output_file_tag];
opt = [opt ' -horizon ' num2str(options_.ms.horizon)];
opt = [opt ' -thin ' num2str(options_.ms.thinning_factor)];
opt = [opt ' -data ' num2str(options_.ms.forecast_data_obs)];

if options_.ms.regimes
    opt = [opt ' -regimes'];
elseif options_.ms.regime
    % regime-1 since regime is 0-indexed in C but 1-indexed in Matlab
    opt = [opt ' -regime ' num2str(options_.ms.regime-1)];
end

if options_.ms.parameter_uncertainty
    options_ = set_ms_simulation_file(options_);
    opt = [opt ' -parameter_uncertainty'];
    opt = [opt ' -shocks_per_parameter ' num2str(options_.ms.shocks_per_parameter)];
else
    opt = [opt ' -shocks_per_parameter ' num2str(options_.ms.shock_draws)];
end

percentiles_size = 0;
if options_.ms.median
    percentiles_size = 1;
    opt = [opt ' -percentiles ' num2str(percentiles_size) ' 0.5'];
else
    percentiles_size = size(options_.ms.percentiles,2);
    opt = [opt ' -percentiles ' num2str(percentiles_size)];
    for i=1:size(options_.ms.percentiles,2)
        opt = [opt ' ' num2str(options_.ms.percentiles(i))];
    end
end

% forecast
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_forecast',err);

% Plot Forecasts
if options_.ms.regimes
    n_chains = length(options_.ms.ms_chain);
    n_regimes=1;
    for i_chain=1:n_chains
        n_regimes = n_regimes*length(options_.ms.ms_chain(i_chain).regime);
    end

    for regime_i=1:n_regimes
        forecast_title = ['Forecast, Regimes ' num2str(regime_i)];
        forecast_data = load([forecastdir filesep 'forecasts_percentiles_regime_' ...
            num2str(regime_i-1) '_' options_.ms.output_file_tag ...
            '.out'], '-ascii');
        forecast_data = reshape_ascii_forecast_data(M_.endo_nbr, ...
            percentiles_size, options_.ms.horizon, forecast_data);
        save([forecastdir filesep 'forecast_regime_' num2str(regime_i-1) '.mat'], ...
            'forecast_data');
        plot_ms_forecast(M_, options_, forecast_data, forecast_title);
    end
else
    if options_.ms.regime
        forecast_data = load([forecastdir filesep 'forecasts_percentiles_regime_' ...
            num2str(options_.ms.regime-1) '_' options_.ms.output_file_tag ...
            '.out'], '-ascii');
        forecast_title = ['Forecast, Regime ' num2str(options_.ms.regime)];
        save_filename = ['forecast_regime_' num2str(options_.ms.regime-1) '.mat'];
    else
        forecast_data = load([forecastdir filesep 'forecasts_percentiles_' ...
            options_.ms.output_file_tag '.out'], '-ascii');
        forecast_title = 'Forecast';
        save_filename = 'forecast.mat';
    end

    forecast_data = reshape_ascii_forecast_data(M_.endo_nbr, ...
        percentiles_size, options_.ms.horizon, forecast_data);
    save([forecastdir filesep save_filename], 'forecast_data');
    plot_ms_forecast(M_, options_, forecast_data, forecast_title);
end
end
