function osr_res = osr1(i_params,i_var,weights)
% Compute the Optimal Simple Rules
% INPUTS
%   i_params                  vector           index of optimizing parameters in M_.params
%   i_var                     vector           variables indices in declaration order
%   weights                   vector           weights in the OSRs
%
% OUTPUTS
%   osr_res:    [structure] results structure containing:
%    - objective_function [scalar double]   value of the objective
%                                               function at the optimum
%    - optim_params       [structure]       parameter values at the optimum 
% 
% Algorithm:
% 
%   Uses Newton-type optimizer csminwel to directly solve quadratic
%   osr-problem
% 
% Copyright (C) 2005-2014 Dynare Team
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

global M_ oo_ options_ it_

klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;


if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
    error ('OSR: Error in model specification: some variables don''t appear as current') ;
end

if M_.maximum_lead == 0
    error ('OSR: Backward or static model: no point in using OSR') ;
end

if any(any(isinf(weights)))
    error ('OSR: At least one of the optim_weights is infinite.') ;
end

if any(isnan(M_.params(i_params)))
    error ('OSR: At least one of the initial parameter values for osr_params is NaN') ;
end

exe =zeros(M_.exo_nbr,1);

oo_.dr = set_state_space(oo_.dr,M_,options_);


np = size(i_params,1);
t0 = M_.params(i_params);


inv_order_var = oo_.dr.inv_order_var;

H0 = 1e-4*eye(np);
crit=options_.osr.tolf;
nit=options_.osr.maxit;

%% do initial checks
[loss,vx,info,exit_flag]=osr_obj(t0,i_params,inv_order_var(i_var),weights(i_var,i_var));
if info~=0
   print_info(info, options_.noprint, options_);
else
   fprintf('\nOSR: Initial value of the objective function: %g \n\n',loss);
end
if isinf(loss)
   fprintf('\nOSR: The initial value of the objective function is infinite.\n');
   fprintf('\nOSR: Check whether the unconditional variance of a target variable is infinite\n');
   fprintf('\nOSR: due to the presence of a unit root.\n');
   error('OSR: Initial likelihood is infinite')
end

%%do actual optimization
[f,p]=csminwel1('osr_obj',t0,H0,[],crit,nit,options_.gradient_method,options_.gradient_epsilon,i_params,...
                inv_order_var(i_var),weights(i_var,i_var));
osr_res.objective_function = f;
for i=1:length(i_params)
    osr_res.optim_params.(deblank(M_.param_names(i_params(i),:))) = p(i);
end

%  options = optimset('fminunc');
%  options = optimset('display','iter');
%  [p,f]=fminunc(@osr_obj,t0,options,i_params,...
%               inv_order_var(i_var),weights(i_var,i_var));

skipline()
disp('OPTIMAL VALUE OF THE PARAMETERS:')
skipline()
for i=1:np
    disp(sprintf('%16s %16.6g\n',M_.param_names(i_params(i),:),p(i)))
end
disp(sprintf('Objective function : %16.6g\n',f));
skipline()
[oo_.dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);