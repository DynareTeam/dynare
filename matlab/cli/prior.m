function varargout = prior(varargin)

% Computes various prior statistics and display them in the command window.
%
% INPUTS
%   'table', 'moments', 'optimize', 'simulate', 'plot'
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2015 Dynare Team
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

if isempty(varargin) || ( isequal(length(varargin), 1) && isequal(varargin{1},'help'))
    skipline()
    disp('Possible options are:')
    disp(' + table       Prints a table describing the priors.')
    disp(' + moments     Computes and displays moments of the endogenous variables at the prior mode.')
    disp(' + optimize    Optimizes the prior density (starting from a random initial guess).')
    disp(' + simulate    Computes the effective prior mass (using a Monte-Carlo).')
    disp(' + plot        Plots the marginal prior densities.')
    skipline()
    return
end

global options_ M_ estim_params_ bayestopt_ oo_ 

donesomething = false;

% Temporarly change qz_criterium option value
changed_qz_criterium_flag  = 0;
if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-9;
    changed_qz_criterium_flag  = 1;
end

M_.dname = M_.fname;

% Temporarly set options_.order equal to one
order = options_.order;
options_.order = 1;

[xparam1,estim_params_,bayestopt_,lb,ub,M_] = set_prior(estim_params_,M_,options_);

if ismember('plot', varargin)
    plot_priors(bayestopt_,M_,estim_params_,options_)
    donesomething = true;
end

if ismember('table', varargin)
    print_table_prior(lb, ub, options_, M_, bayestopt_, estim_params_);
    donesomething = true;
end

if ismember('simulate', varargin) % Prior simulations (BK).
    results = prior_sampler(0,M_,bayestopt_,options_,oo_,estim_params_);
    % Display prior mass info
    skipline(2)
    disp(['Prior mass = ' num2str(results.prior.mass)])
    disp(['BK indeterminacy share                = ' num2str(results.bk.indeterminacy_share)])
    disp(['BK unstability share                  = ' num2str(results.bk.unstability_share)])
    disp(['BK singularity share                  = ' num2str(results.bk.singularity_share)])
    disp(['Complex jacobian share                = ' num2str(results.jacobian.problem_share)])
    disp(['mjdgges crash share                   = ' num2str(results.dll.problem_share)])
    disp(['Steady state problem share            = ' num2str(results.ss.problem_share)])
    disp(['Complex steady state  share           = ' num2str(results.ss.complex_share)])
    disp(['Analytical steady state problem share = ' num2str(results.ass.problem_share)])
    skipline(2)
    donesomething = true;
end

if ismember('optimize', varargin) % Prior optimization.
    optimize_prior(options_, M_, oo_, bayestopt_, estim_params_);
    donesomething = true;
end

if ismember('moments', varargin) % Prior simulations (2nd order moments).
    % Set estimated parameters to the prior mode...
    xparam1 = bayestopt_.p5;
    % ... Except for uniform priors!
    k = find(isnan(xparam1));
    xparam1(k) = bayestopt_.p5(k);
    % Update vector of parameters and covariance matrices
    M_ = set_all_parameters(xparam1,estim_params_,M_);
    % Check model.
    check_model(M_);
    % Compute state space representation of the model.
    oo_.dr=set_state_space(oo_.dr, M_, options_);
    % Solve model
    [dr,info, M_ ,options_ , oo_] = resol(0, M_ , options_ ,oo_);
    % Compute and display second order moments
    oo_=disp_th_moments(oo_.dr,[],M_,options_,oo_);
    skipline(2)
    donesomething = true;
end

if changed_qz_criterium_flag
    options_.qz_criterium = [];
end

options_.order = order;

if ~donesomething
    error('prior: Unexpected arguments!')
end

function format_string = build_format_string(PriorStandardDeviation,LowerBound,UpperBound)
format_string = ['%s & %s & %6.4f &'];
if ~isnumeric(PriorStandardDeviation)
    format_string = [ format_string , ' %s &'];
else
    format_string = [ format_string , ' %6.4f &'];
end
if ~isnumeric(LowerBound)
    format_string = [ format_string , ' %s &'];
else
    format_string = [ format_string , ' %6.4f &'];
end
if ~isnumeric(UpperBound)
    format_string = [ format_string , ' %s &'];
else
    format_string = [ format_string , ' %6.4f &'];
end
format_string = [ format_string , ' %6.4f & %6.4f \\\\ \n'];