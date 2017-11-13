function oo_ = initial_condition_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% function oo_ = initial_condition_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% Computes initial condition contribution to a simulated trajectory. The field set is
% oo_.initval_decomposition. It is a n_var by n_var+2 by nperiods array. The
% first n_var columns store the respective endogenous initval contribution, column n+1
% stores the role of the shocks, while column n+2 stores the
% value of the smoothed variables.  Variables are stored
% in the order of declaration, i.e. M_.endo_names.
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%    bayestopt_:  [structure]  describing the priors
%    estim_params_: [structure] characterizing parameters to be estimated
%
% OUTPUTS
%    oo_:         [structure]  Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017 Dynare Team
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

options_.plot_shock_decomp.detail_plot = options_.initial_condition_decomp.detail_plot;
options_.plot_shock_decomp.steadystate = options_.initial_condition_decomp.steadystate;
options_.plot_shock_decomp.write_xls = options_.initial_condition_decomp.write_xls;
options_.plot_shock_decomp.type = options_.initial_condition_decomp.type;
options_.plot_shock_decomp.plot_init_date = options_.initial_condition_decomp.plot_init_date;
options_.plot_shock_decomp.plot_end_date = options_.initial_condition_decomp.plot_end_date;

% indices of endogenous variables
if size(varlist,1) == 0
    varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
end

[i_var,nvar,index_uniques] = varlist_indices(varlist,M_.endo_names);
varlist=varlist(index_uniques,:);

% number of variables
endo_nbr = M_.endo_nbr;

% parameter set
parameter_set = options_.parameter_set;
if isempty(parameter_set)
    if isfield(oo_,'posterior_mean')
        parameter_set = 'posterior_mean';
    elseif isfield(oo_,'mle_mode')
        parameter_set = 'mle_mode';
    elseif isfield(oo_,'posterior')
        parameter_set = 'posterior_mode';
    else
        error(['shock_decomposition: option parameter_set is not specified ' ...
               'and posterior mode is not available'])
    end
end

if ~isfield(oo_,'initval_decomposition')
    options_.selected_variables_only = 0; %make sure all variables are stored
    options_.plot_priors=0;
    [oo,M,junk1,junk2,Smoothed_Variables_deviation_from_mean] = evaluate_smoother(parameter_set,varlist,M_,oo_,options_,bayestopt_,estim_params_);

    % reduced form
    dr = oo.dr;

    % data reordering
    order_var = dr.order_var;
    inv_order_var = dr.inv_order_var;


    % coefficients
    A = dr.ghx;
    B = dr.ghu;

    % initialization
    gend = size(oo.SmoothedShocks.(deblank(M_.exo_names(1,:))),1); %+options_.forecast;
    z = zeros(endo_nbr,endo_nbr+2,gend);
    z(:,end,:) = Smoothed_Variables_deviation_from_mean;

    for i=1:endo_nbr
        z(i,i,1) = Smoothed_Variables_deviation_from_mean(i,1);
    end

    maximum_lag = M_.maximum_lag;

    k2 = dr.kstate(find(dr.kstate(:,2) <= maximum_lag+1),[1 2]);
    i_state = order_var(k2(:,1))+(min(i,maximum_lag)+1-k2(:,2))*M_.endo_nbr;
    for i=1:gend
        if i > 1 && i <= maximum_lag+1
            lags = min(i-1,maximum_lag):-1:1;
        end

        if i > 1
            tempx = permute(z(:,1:endo_nbr,lags),[1 3 2]);
            m = min(i-1,maximum_lag);
            tempx = [reshape(tempx,endo_nbr*m,endo_nbr); zeros(endo_nbr*(maximum_lag-i+1),endo_nbr)];
            z(:,1:endo_nbr,i) = A(inv_order_var,:)*tempx(i_state,:);
            lags = lags+1;
        end
        z(:,endo_nbr+1,i) = z(:,endo_nbr+2,i) - sum(z(:,1:endo_nbr,i),2);

    end


    oo_.initval_decomposition = z;
end
% if ~options_.no_graph.shock_decomposition
oo=oo_;
oo.shock_decomposition = oo_.initval_decomposition;
M_.exo_names = M_.endo_names;
M_.exo_nbr = M_.endo_nbr;
options_.plot_shock_decomp.realtime=0;
options_.plot_shock_decomp.screen_shocks=1;
options_.plot_shock_decomp.use_shock_groups = '';
options_.plot_shock_decomp.fig_name='initval';
plot_shock_decomposition(M_,oo,options_,varlist);
% end