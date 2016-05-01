function [flag,endo_simul,err,y] = solve_stochastic_perfect_foresight_model_1(endo_simul,exo_simul,Options,pfm,order,varargin)

% Copyright (C) 2012-2013 Dynare Team
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

if nargin < 6
    homotopy_parameter = 1;
else
    homotopy_parameter = varargin{1};
end
flag = 0;
err = 0;
stop = 0;

EpOptions = Options.ep;

params = pfm.params;
steady_state = pfm.steady_state;
ny = pfm.ny;
periods = pfm.periods;
dynamic_model = pfm.dynamic_model;
lead_lag_incidence = pfm.lead_lag_incidence;
nyp = pfm.nyp;
nyf = pfm.nyf;
i_cols_1 = pfm.i_cols_1;
i_cols_A1 = pfm.i_cols_A1;
i_cols_j = pfm.i_cols_j;
i_cols_T = nonzeros(lead_lag_incidence(1:2,:)');
hybrid_order = pfm.hybrid_order;
dr = pfm.dr;
nodes = pfm.nodes;
weights = pfm.weights;
nnodes = pfm.nnodes;

maxit = pfm.maxit_;
tolerance = pfm.tolerance;
verbose = pfm.verbose;

number_of_shocks = size(exo_simul,2);

% make sure that there is a node equal to zero
% and permute nodes and weights to have zero first
k = find(sum(abs(nodes),2) < 1e-12);
if ~isempty(k)
    nodes = [nodes(k,:); nodes(1:k-1,:); nodes(k+1:end,:)];
    weights = [weights(k); weights(1:k-1); weights(k+1:end)];
else
    error('there is no nodes equal to zero')
end
if hybrid_order > 0
    if hybrid_order == 2
        h_correction = 0.5*dr.ghs2(dr.inv_order_var);
    end
else
    h_correction = 0;
end

if verbose
    disp ([' -----------------------------------------------------']);
    disp (['MODEL SIMULATION :']);
    fprintf('\n');
end

% Each column of Y represents a different world
% The upper right cells are unused
% The first row block is ny x 1
% The second row block is ny x nnodes
% The third row block is ny x nnodes^2
% and so on until size ny x nnodes^order
world_nbr = pfm.world_nbr;
Y = endo_simul(:,2:end-1);
Y = repmat(Y,1,world_nbr); 
pfm.y0 = endo_simul(:,1);

% The columns of A map the elements of Y such that
% each block of Y with ny rows are unfolded column wise
% number of blocks
block_nbr = pfm.block_nbr;
% dimension of the problem
dimension = ny*block_nbr;
pfm.dimension = dimension;
if order == 0
    i_upd_r = (1:ny*periods)';
    i_upd_y = i_upd_r + ny;
else
    i_upd_r = zeros(dimension,1);
    i_upd_y = i_upd_r;
    i_upd_r(1:ny) = (1:ny);
    i_upd_y(1:ny) = ny+(1:ny);
    i1 = ny+1;
    i2 = 2*ny;
    n1 = ny+1;
    n2 = 2*ny;
    for i=2:periods
        k = n1:n2;
        for j=1:(1+(nnodes-1)*min(i-1,order))
            i_upd_r(i1:i2) = k+(j-1)*ny*periods;
            i_upd_y(i1:i2) = k+ny+(j-1)*ny*(periods+2);
            i1 = i2+1;
            i2 = i2+ny;
        end
        n1 = n2+1;
        n2 = n2+ny;
    end
end
icA = [find(lead_lag_incidence(1,:)) find(lead_lag_incidence(2,:))+world_nbr*ny ...
       find(lead_lag_incidence(3,:))+2*world_nbr*ny]';
h1 = clock;
pfm.order = order;
pfm.world_nbr = world_nbr;
pfm.nodes = nodes;
pfm.nnodes = nnodes;
pfm.weights = weights;
pfm.h_correction = h_correction;
pfm.i_rows = 1:ny;
i_cols = find(lead_lag_incidence');
pfm.i_cols = i_cols;
pfm.nyp = nyp;
pfm.nyf = nyf;
pfm.hybrid_order = hybrid_order;
pfm.i_cols_1 = i_cols_1;
pfm.i_cols_h = i_cols_j;
pfm.icA = icA;
pfm.i_cols_T = i_cols_T;
pfm.i_upd_r = i_upd_r;
pfm.i_upd_y = i_upd_y;

y = repmat(steady_state,block_nbr,1);
Options.solve_algo = Options.ep.solve_algo;
Options.steady.maxit = Options.ep.maxit;
[y,info] = dynare_solve(@ep_problem_2,y,Options,exo_simul,pfm);
if info
    flag = 1;
    err = info;
end
endo_simul(:,2) = y(1:ny);