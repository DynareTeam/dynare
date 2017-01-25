function osr_res = osr(var_list,params,i_var,W)
% osr_res = osr(var_list,params,i_var,W)
%   Wrapper function computing the solution to the optimal simple
%   rule-problem; calls osr1 for actual computation
% INPUTS
%   var_list    [character array]           list of endogenous variables specified
%                                           after osr1-command (deprecated and not used anymore)
%   params      [character array]           list of parameter to be chosen in
%                                           optimal simple rule
%   i_var       [n_osr_vars by 1 double]    indices of osr-variable in
%                                           specified in optim_weights in declaration order
%   W           [M_.endo_nbr by M_.endo_nbr sparse matrix] Weighting matrix for variance of endogenous variables             
%    
% OUTPUTS
%   osr_res:    [structure] results structure containing:
%    - objective_function [scalar double]   value of the objective
%                                               function at the optimum
%    - optim_params       [structure]       parameter values at the optimum 
% 
%       
% SPECIAL REQUIREMENTS
%   none.
%  
% Copyright (C) 2001-2016 Dynare Team
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

options_.order = 1;

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

oo_=make_ex_(M_,options_,oo_);

np = size(params,1);
i_params = zeros(np,1);
for i=1:np
    str = deblank(params(i,:));
    i_params(i) = strmatch(str{:}, M_.param_names, 'exact');
end

if ~options_.noprint
    skipline()
    disp('OPTIMAL SIMPLE RULE')
    skipline()
end
osr_res = osr1(i_params,i_var,W);

stoch_simul(var_list);