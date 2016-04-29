function [ts, DynareResults] = extended_path(initialconditions, samplesize, exogenousvariables, DynareOptions, DynareModel, DynareResults)

% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS 
%  o initialconditions      [double]    m*1 array, where m is the number of endogenous variables in the model.
%  o samplesize             [integer]   scalar, size of the sample to be simulated.
%  o exogenousvariables     [double]    T*n array, values for the structural innovations.
%  o DynareOptions          [struct]    options_
%  o DynareModel            [struct]    M_
%  o DynareResults          [struct]    oo_
%
% OUTPUTS 
%  o ts                     [dseries]   m*samplesize array, the simulations.
%  o results                [cell]
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2009-2016 Dynare Team
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

[initialconditions, innovations, pfm, ep, verbosity, DynareOptions, DynareResults] = ...
    extended_path_initialization(initialconditions, samplesize, exogenousvariables, DynareOptions, DynareModel, DynareResults);

[shocks, spfm_exo_simul, innovations, DynareResults] = extended_path_shocks(innovations, ep, exogenousvariables, samplesize, DynareResults);

% Initialize the matrix for the paths of the endogenous variables.
endogenous_variables_paths = NaN(DynareModel.endo_nbr,samplesize+1);
endogenous_variables_paths(:,1) = initialconditions;

% Set waitbar (graphic or text  mode)
hh = dyn_waitbar(0,'Please wait. Extended Path simulations...');
set(hh,'Name','EP simulations.');

% Initialize while-loop index.
t = 1;

% Main loop.
while (t <= samplesize)
    if ~mod(t,10)
        dyn_waitbar(t/samplesize,hh,'Please wait. Extended Path simulations...');
    end
    % Set period index.
    t = t+1;
    spfm_exo_simul(2,:) = shocks(t-1,:);
    if t>2
        % Set initial guess for the solver (using the solution of the
        % previous period problem).
        initialguess = [endogenousvariablespaths(:, 2:end), DynareResults.steady_state];
    else
        initialguess = [];
    end
    [endogenous_variables_paths(:,t), info_convergence, endogenousvariablespaths] = extended_path_core(ep.periods, DynareModel.endo_nbr, DynareModel.exo_nbr, innovations.positive_var_indx, ...
                                                      spfm_exo_simul, ep.init, endogenous_variables_paths(:,t-1), ...
                                                      DynareResults.steady_state, ...
                                                      ep.verbosity, ep.use_bytecode, ep.stochastic.order, ...
                                                      DynareModel, pfm, ep.stochastic.algo, ep.solve_algo, ep.stack_solve_algo, ...
                                                      DynareOptions.lmmcp, ...
                                                      DynareOptions, ...
                                                      DynareResults, initialguess);
    if ~info_convergence
        msg = sprintf('No convergence of the (stochastic) perfect foresight solver (in period %s)!', int2str(t));
        warning(msg)
        break
    end
end% (while) loop over t

% Close waitbar.
dyn_waitbar_close(hh);

% Set the initial period.
if isnan(DynareOptions.initial_period)
    initial_period = dates(1,1);
else
    initial_period = DynareOptions.initial_period;
end

% Return the simulated time series.
if any(isnan(endogenous_variables_paths(:)))
    sl = find(~isnan(endogenous_variables_paths));
    nn = size(endogenous_variables_paths, 1);
    endogenous_variables_paths = reshape(endogenous_variables_paths(sl), nn, length(sl)/nn); 
end
ts = dseries(transpose(endogenous_variables_paths), initial_period, cellstr(DynareModel.endo_names));

DynareResults.endo_simul = transpose(ts.data);
assignin('base', 'Simulated_time_series', ts);

if ~nargout || nargout<2
    assignin('base', 'oo_', DynareResults);
end