function [options_, oo_]=ms_irf(varlist, M_, options_, oo_)
% function [options_, oo_]=ms_irf(varlist, M_, options_, oo_)
% Markov-switching SBVAR: Impulse Response Function
%
% INPUTS
%    varlist:     (chararray) list of selected endogenous variables
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

disp('MS-SBVAR Impulse Response Function');
options_ = set_file_tags(options_);
[options_, oo_] = set_ms_estimation_file(options_.ms.file_tag, options_, oo_);
clean_ms_irf_files(options_.ms.output_file_tag);
irfdir = [options_.ms.output_file_tag filesep 'IRF'];
create_dir(irfdir);

% setup command line options
opt = ['-ir -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -do ' irfdir];
opt = [opt ' -ft ' options_.ms.file_tag];
opt = [opt ' -fto ' options_.ms.output_file_tag];
opt = [opt ' -horizon ' num2str(options_.ms.horizon)];
opt = [opt ' -thin ' num2str(options_.ms.thinning_factor)];

if options_.ms.regimes
    opt = [opt ' -regimes'];
elseif options_.ms.regime
    % regime-1 since regime is 0-indexed in C but 1-indexed in Matlab
    opt = [opt ' -regime ' num2str(options_.ms.regime-1)];
elseif options_.ms.filtered_probabilities
    opt = [opt ' -filtered'];
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

% irf
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_irf',err);

% Plot IRFs
if options_.ms.regimes
    n_chains = length(options_.ms.ms_chain);
    n_regimes=1;
    for i_chain=1:n_chains
        n_regimes = n_regimes*length(options_.ms.ms_chain(i_chain).regime);
    end

    for regime_i=1:n_regimes
        irf_title = ['Impulse Responses, Regime ' num2str(regime_i)];
        irf_data = load([irfdir filesep 'ir_percentiles_regime_' ...
            num2str(regime_i-1) '_' options_.ms.output_file_tag ...
            '.out'], '-ascii');
        irf_data = reshape_ascii_irf_data(M_.endo_nbr, percentiles_size, ...
            options_.ms.horizon, irf_data);
        save([irfdir filesep 'irf_regime_' num2str(regime_i-1) '.mat'], 'irf_data');
        plot_ms_irf(M_, options_, irf_data, irf_title, varlist);
    end
else
    if options_.ms.regime
        irf_data = load([irfdir filesep 'ir_percentiles_regime_' ...
            num2str(options_.ms.regime-1) '_' options_.ms.output_file_tag ...
            '.out'], '-ascii');
        irf_title = ['Impulse Response, Regime ' num2str(options_.ms.regime)];
        save_filename = ['irf_regime_' num2str(options_.ms.regime-1) '.mat'];
    elseif options_.ms.filtered_probabilities
        irf_data = load([irfdir filesep 'ir_percentiles_filtered_' ...
            options_.ms.output_file_tag '.out'], '-ascii');
        irf_title = 'Impulse Response Filtered';
        save_filename = 'irf.mat';
    else
        irf_data = load([irfdir filesep 'ir_percentiles_ergodic_' ...
            options_.ms.output_file_tag '.out'], '-ascii');
        irf_title = 'Impulse Response Ergodic';
        save_filename = 'irf.mat';
    end

    irf_data = reshape_ascii_irf_data(M_.endo_nbr, percentiles_size, ...
        options_.ms.horizon, irf_data);
    save([irfdir filesep save_filename], 'irf_data');
    plot_ms_irf(M_, options_, irf_data, irf_title, varlist);
end
end
