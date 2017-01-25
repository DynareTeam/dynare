function [options_, oo_]=ms_variance_decomposition(M_, options_, oo_)
% function [options_, oo_]=ms_variance_decomposition(M_, options_, oo_)
% Markov-switching SBVAR: Variance Decomposition
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

disp('MS-SBVAR Variance Decomposition');
options_ = set_file_tags(options_);
[options_, oo_] = set_ms_estimation_file(options_.ms.file_tag, options_, oo_);
clean_ms_variance_decomposition_files(options_.ms.output_file_tag);
vddir = [options_.ms.output_file_tag filesep 'Variance_Decomposition'];
create_dir(vddir);

% setup command line options
opt = ['-variance_decomposition -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -do ' vddir];
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

percentiles_size = 1;
outfile = [vddir filesep 'var_decomp_mean_'];
if options_.ms.error_bands
    % error_bands / percentiles used differently by
    % Dan's variance decomposition code
    % no_error_bands => mean is computed
    percentiles_size = size(options_.ms.percentiles,2);
    opt = [opt ' -percentiles ' num2str(percentiles_size)];
    for i=1:size(options_.ms.percentiles,2)
        opt = [opt ' ' num2str(options_.ms.percentiles(i))];
    end
    outfile = [vddir filesep 'var_decomp_percentiles_'];
end

% variance_decomposition
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_variance_decomposition',err);

if options_.ms.regime || options_.ms.regimes
    outfile = [outfile 'regime_'];
    if options_.ms.regime
        outfile = [outfile num2str(options_.ms.regime-1) ...
            '_' options_.ms.output_file_tag '.out'];
    end
elseif options_.ms.filtered_probabilities
    outfile = [outfile 'filtered_' options_.ms.output_file_tag '.out'];
else
    outfile = [outfile 'ergodic_' options_.ms.output_file_tag '.out'];
end

% Create plots
if options_.ms.regimes
    n_chains = length(options_.ms.ms_chain);
    n_regimes=1;
    for i_chain=1:n_chains
        n_regimes = n_regimes*length(options_.ms.ms_chain(i_chain).regime);
    end
    for regime_i=1:n_regimes
        vd_title = ['Variance Decomposition, Regime ' num2str(regime_i)];
        vd_data = load([outfile num2str(regime_i-1) '_' ...
            options_.ms.output_file_tag '.out'], '-ascii');
        vd_data = reshape_ascii_variance_decomposition_data( ...
            M_.endo_nbr, percentiles_size, options_.ms.horizon, vd_data);
        save([vddir filesep 'variance_decomposition_regime_' num2str(regime_i-1) '.mat'], 'vd_data');
        plot_ms_variance_decomposition(M_, options_, vd_data, vd_title);
    end
else
    if options_.ms.regime
        vd_title = ['Variance Decomposition, Regime ' num2str(options_.ms.regime)];
        save_filename = ['variance_decomposition_regime_' num2str(options_.ms.regime-1) '.mat'];
    else
        save_filename = 'variance_decomposition.mat';
        if options_.ms.filtered_probabilities
            vd_title = 'Variance Decomposition Filtered';
        else
            vd_title = 'Variance Decomposition Ergodic';
        end
    end
    vd_data = load(outfile, '-ascii');
    vd_data = reshape_ascii_variance_decomposition_data( ...
        M_.endo_nbr, percentiles_size, options_.ms.horizon, vd_data);
    save([vddir filesep save_filename], 'vd_data');
    plot_ms_variance_decomposition(M_, options_, vd_data, vd_title);
end
end
