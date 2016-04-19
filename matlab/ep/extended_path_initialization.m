function [initial_conditions, innovations, pfm, ep, verbosity, DynareOptions, DynareResults] = extended_path_initialization(initial_conditions, sample_size, exogenousvariables, DynareOptions, DynareModel, DynareResults)

% Initialization of the extended path routines.
%
% INPUTS 
%  o initial_conditions     [double]    m*1 array, where m is the number of endogenous variables in the model.
%  o sample_size            [integer]   scalar, size of the sample to be simulated.
%  o exogenousvariables     [double]    T*n array, values for the structural innovations.
%  o DynareOptions          [struct]    options_ 
%  o DynareModel            [struct]    M_
%  o DynareResults          [struct]    oo_
%
% OUTPUTS 
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

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

ep  = DynareOptions.ep;

% Set verbosity levels.
DynareOptions.verbosity = ep.verbosity;
verbosity = ep.verbosity+ep.debug;

% Set maximum number of iterations for the deterministic solver.
DynareOptions.simul.maxit = ep.maxit;

% Prepare a structure needed by the matlab implementation of the perfect foresight model solver
pfm = setup_stochastic_perfect_foresight_model_solver(DynareModel, DynareOptions, DynareResults);

% Check that the user did not use varexo_det
if DynareModel.exo_det_nbr~=0
    error('Extended path does not support varexo_det.')
end

% Set default initial conditions.
if isempty(initial_conditions)
    if isempty(DynareModel.endo_histval)
        initial_conditions = DynareResults.steady_state;
    else
        initial_conditions = DynareModel.endo_histval;
    end
end

% Set the number of periods for the (stochastic) perfect foresight model
pfm.periods = ep.periods;

pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);

pfm.block = DynareOptions.block;

% Set the algorithm for the perfect foresight solver
DynareOptions.stack_solve_algo = ep.stack_solve_algo;

% Compute the first order reduced form if needed.
%
% REMARK. It is assumed that the user did run the same mod file with stoch_simul(order=1) and save
% all the globals in a mat file called linear_reduced_form.mat;

dr = struct();
if ep.init
    DynareOptions.order = 1;
    DynareResults.dr=set_state_space(dr,DynareModel,DynareOptions);
    [dr,Info,DynareModel,DynareOptions,DynareResults] = resol(0,DynareModel,DynareOptions,DynareResults);
end

% Do not use a minimal number of perdiods for the perfect foresight solver (with bytecode and blocks)
DynareOptions.minimal_solving_period = DynareOptions.ep.periods;

% Set the covariance matrix of the structural innovations.
if isempty(exogenousvariables)
    innovations = struct();
    innovations.positive_var_indx = find(diag(DynareModel.Sigma_e)>0);
    innovations.effective_number_of_shocks = length(innovations.positive_var_indx);
    innovations.covariance_matrix = DynareModel.Sigma_e(innovations.positive_var_indx,innovations.positive_var_indx);
    innovations.covariance_matrix_upper_cholesky = chol(innovations.covariance_matrix);
else
    innovations = struct();
end

% Set seed.
if ep.set_dynare_seed_to_default
    set_dynare_seed('default');
end

% hybrid correction
pfm.hybrid_order = ep.stochastic.hybrid_order;
if pfm.hybrid_order
    DynareResults.dr = set_state_space(DynareResults.dr, DynareModel, DynareOptions);
    options = DynareOptions;
    options.order = pfm.hybrid_order;
    pfm.dr = resol(0, DynareModel, options, DynareResults);
else
    pfm.dr = [];
end

% number of nonzero derivatives
pfm.nnzA = DynareModel.NNZDerivatives(1);

% setting up integration nodes if order > 0
if ep.stochastic.order > 0
    [nodes,weights,nnodes] = setup_integration_nodes(DynareOptions.ep,pfm);
    pfm.nodes = nodes;
    pfm.weights = weights; 
    pfm.nnodes = nnodes;
    % compute number of blocks
    [block_nbr,pfm.world_nbr] = get_block_world_nbr(ep.stochastic.algo,nnodes,ep.stochastic.order,ep.periods);
else
    block_nbr = ep.periods;
end

% set boundaries if mcp
[lb,ub,pfm.eq_index] = get_complementarity_conditions(DynareModel, DynareOptions.ramsey_policy);
DynareOptions.lmmcp.lb = repmat(lb,block_nbr,1);
DynareOptions.lmmcp.ub = repmat(ub,block_nbr,1);
pfm.block_nbr = block_nbr;

