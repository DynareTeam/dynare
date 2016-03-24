function [ts,results] = extended_path(initial_conditions,sample_size)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS
%  o initial_conditions     [double]    m*nlags array, where m is the number of endogenous variables in the model and
%                                       nlags is the maximum number of lags.
%  o sample_size            [integer]   scalar, size of the sample to be simulated.
%
% OUTPUTS
%  o time_series            [double]    m*sample_size array, the simulations.
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
    global M_ options_ oo_

ep  = options_.ep;
options_.verbosity = ep.verbosity;
verbosity = ep.verbosity+ep.debug;

% Set maximum number of iterations for the deterministic solver.
options_.simul.maxit = ep.maxit;

% Prepare a structure needed by the matlab implementation of the perfect foresight model solver
pfm = setup_stochastic_perfect_foresight_model_solver(M_,options_,oo_);

if M_.exo_det_nbr~=0
    error('ep: Extended path does not support varexo_det.')
end

endo_nbr = M_.endo_nbr;
exo_nbr = M_.exo_nbr;
maximum_lag = M_.maximum_lag;
maximum_lead = M_.maximum_lead;
epreplic_nbr = ep.replic_nbr;
steady_state = oo_.steady_state;
dynatol = options_.dynatol;

% Set default initial conditions.
if isempty(initial_conditions)
    if isempty(M_.endo_histval)
        initial_conditions = steady_state;
    else
        initial_conditions = M_.endo_histval;
    end
end


% Set the number of periods for the perfect foresight model
periods = ep.periods;
pfm.periods = periods;
pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);
pfm.block = options_.block;

% keep a copy of pfm.i_upd
i_upd = pfm.i_upd;

% Set the algorithm for the perfect foresight solver
options_.stack_solve_algo = ep.stack_solve_algo;

% Set check_stability flag
do_not_check_stability_flag = ~ep.check_stability;

% Compute the first order reduced form if needed.
%
% REMARK. It is assumed that the user did run the same mod file with stoch_simul(order=1) and save
% all the globals in a mat file called linear_reduced_form.mat;

dr = struct();
if ep.init
    options_.order = 1;
    oo_.dr=set_state_space(dr,M_,options_);
    [dr,Info,M_,options_,oo_] = resol(0,M_,options_,oo_);
end

% Do not use a minimal number of perdiods for the perfect foresight solver (with bytecode and blocks)
options_.minimal_solving_period = 100;%options_.ep.periods;

% Initialize the output array.
time_series = zeros(M_.endo_nbr,sample_size);

% Set the covariance matrix of the structural innovations.
variances = diag(M_.Sigma_e);
positive_var_indx = find(variances>0);
effective_number_of_shocks = length(positive_var_indx);
stdd = sqrt(variances(positive_var_indx));
covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx);
covariance_matrix_upper_cholesky = chol(covariance_matrix);

% (re)Set exo_nbr
%exo_nbr = effective_number_of_shocks;

% Set seed.
if ep.set_dynare_seed_to_default
    set_dynare_seed('default');
end

% Set bytecode flag
bytecode_flag = ep.use_bytecode;
% Set number of replications
replic_nbr = ep.replic_nbr;

% Simulate shocks.
switch ep.innovation_distribution
  case 'gaussian'
    shocks = transpose(transpose(covariance_matrix_upper_cholesky)* ...
                       randn(effective_number_of_shocks,sample_size* ...
                             replic_nbr));
    shocks(:,positive_var_indx) = shocks;
  case 'calibrated'
    replic_nbr = 1;
    shocks = zeros(sample_size,M_.exo_nbr);
    for i = 1:length(M_.unanticipated_det_shocks)
        k = M_.unanticipated_det_shocks(i).periods;
        ivar = M_.unanticipated_det_shocks(i).exo_id;
        v = M_.unanticipated_det_shocks(i).value;
        if ~M_.unanticipated_det_shocks(i).multiplicative
            shocks(k,ivar) = v;
        else
            socks(k,ivar) = shocks(k,ivar) * v;
        end
    end
  otherwise
    error(['extended_path:: ' ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
end

% Set waitbar (graphic or text  mode)
hh = dyn_waitbar(0,'Please wait. Extended Path simulations...');
set(hh,'Name','EP simulations.');

% hybrid correction
pfm.hybrid_order = ep.stochastic.hybrid_order;
if pfm.hybrid_order
    oo_.dr = set_state_space(oo_.dr,M_,options_);
    options = options_;
    options.order = pfm.hybrid_order;
    pfm.dr = resol(0,M_,options,oo_);
else
    pfm.dr = [];
end

% number of nonzero derivatives
pfm.nnzA = M_.NNZDerivatives(1);

% setting up integration nodes if order > 0
if ep.stochastic.order > 0
    [nodes,weights,nnodes] = setup_integration_nodes(options_.ep,pfm);
    pfm.nodes = nodes;
    pfm.weights = weights; 
    pfm.nnodes = nnodes;

    % compute number of blocks
    [block_nbr,pfm.world_nbr] = get_block_world_nbr(ep.stochastic.algo,nnodes,ep.stochastic.order,ep.periods);
else
    block_nbr = ep.periods;
end


% set boundaries if mcp
[lb,ub,pfm.eq_index] = get_complementarity_conditions(M_, options_.ramsey_policy);
options_.lmmcp.lb = repmat(lb,block_nbr,1);
options_.lmmcp.ub = repmat(ub,block_nbr,1);
pfm.block_nbr = block_nbr;

% storage for failed draws
oo_.ep.failures.periods = [];
oo_.ep.failures.previous_period = cell(0);
oo_.ep.failures.shocks = cell(0);

oo_.exo_simul = shocks;

% Initializes some variables.
t  = 1;
tsimul = 1;
for k = 1:replic_nbr
    results{k} = zeros(endo_nbr,sample_size+1);
    results{k}(:,1) = initial_conditions;
end
%make_ex_;
exo_simul_ = zeros(maximum_lag+sample_size+maximum_lead,exo_nbr);
exo_simul_(1:size(oo_.exo_simul,1),1:size(oo_.exo_simul,2)) = oo_.exo_simul;
% Main loop.
while (t <= sample_size)
    if ~mod(t,10)
        dyn_waitbar(t/sample_size,hh,'Please wait. Extended Path simulations...');
    end
    % Set period index.
    t = t+1;
    info_convergence=1;
    if replic_nbr > 1 && ep.parallel_1
        parfor k = 1:replic_nbr
            exo_simul = repmat(oo_.exo_steady_state',periods+2,1);
            %            exo_simul(1:sample_size+3-t,:) = exo_simul_(t:end,:);
            exo_simul(2,:) = exo_simul_(M_.maximum_lag+t,:) + ...
                shocks((t-2)*replic_nbr+k,:);
            initial_conditions = results{k}(:,t-1);
            [results{k}(:,t), tmp_info_convergence] = extended_path_core(ep.periods,endo_nbr,exo_nbr,positive_var_indx, ...
                                                              exo_simul,ep.init,initial_conditions,...
                                                              maximum_lag,maximum_lead,steady_state, ...
                                                              ep.verbosity,bytecode_flag,ep.stochastic.order,...
                                                              M_.params,pfm,ep.stochastic.algo,ep.solve_algo,ep.stack_solve_algo,...
                                                              options_.lmmcp,options_,oo_);
            info_convergence = min(tmp_info_convergence, info_convergence);                                              
        end
    else
        for k = 1:replic_nbr
            exo_simul = repmat(oo_.exo_steady_state',periods+maximum_lag+ ...
                            maximum_lead,1);
            %            exo_simul(1:sample_size+maximum_lag+maximum_lead-t+1,:) = ...
            %                exo_simul_(t:end,:);
            exo_simul(maximum_lag+1,:) = ...
                exo_simul_(maximum_lag+t,:) + shocks((t-2)*replic_nbr+k,:);
            initial_conditions = results{k}(:,t-1);
            [results{k}(:,t), info_convergence] = extended_path_core(ep.periods,endo_nbr,exo_nbr,positive_var_indx, ...
                                                              exo_simul,ep.init,initial_conditions,...
                                                              maximum_lag,maximum_lead,steady_state, ...
                                                              ep.verbosity,bytecode_flag,ep.stochastic.order,...
                                                              M_,pfm,ep.stochastic.algo,ep.solve_algo,ep.stack_solve_algo,...
                                                              options_.lmmcp,options_,oo_);
        end
    end
    if verbosity
        if info_convergence
            disp(['Time: ' int2str(t)  '. Convergence of the perfect foresight model solver!'])
        else
            disp(['Time: ' int2str(t)  '. No convergence of the perfect foresight model solver!'])
        end
    end
end% (while) loop over t

dyn_waitbar_close(hh);

if isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    initial_period = options_.initial_period;
end
if nargout
    if ~isnan(results{1})
        ts = dseries(transpose([results{1}]), ...
                     initial_period,cellstr(M_.endo_names));
    else
        ts = NaN;
    end
else
    if ~isnan(results{1})
        oo_.endo_simul = results{1};
        ts = dseries(transpose(results{1}),initial_period, ...
                     cellstr(M_.endo_names));
    else
        oo_.endo_simul = NaN;
        ts = NaN;
    end
end

 assignin('base', 'Simulated_time_series', ts);