function oo_ = display_conditional_variance_decomposition(Steps, SubsetOfVariables, dr,M_,options_,oo_)
% This function computes the conditional variance decomposition of a given state space model
% for a subset of endogenous variables.
% 
% INPUTS 
%   StateSpaceModel     [structure]   Specification of the state space model.
%   Steps               [integer]     1*h vector of dates.
%   SubsetOfVariables   [integer]     1*q vector of indices.
%    
% OUTPUTS 
%   PackedConditionalVarianceDecomposition  [double] n(n+1)/2*p matrix, where p is the number of state innovations and
%                                                    n is equal to length(SubsetOfVariables).    
%
% SPECIAL REQUIREMENTS
%
% [1] The covariance matrix of the state innovations needs to be diagonal.
% [2] In this version, absence of measurement errors is assumed...

% Copyright (C) 2010-2012 Dynare Team
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

endo_nbr = M_.endo_nbr;
exo_nbr = M_.exo_nbr;
StateSpaceModel.number_of_state_equations = M_.endo_nbr;
StateSpaceModel.number_of_state_innovations = exo_nbr;
StateSpaceModel.sigma_e_is_diagonal = M_.sigma_e_is_diagonal;

iv = (1:endo_nbr)';
ic = M_.nstatic+(1:M_.nspred)';

[StateSpaceModel.transition_matrix,StateSpaceModel.impulse_matrix] = kalman_transition_matrix(dr,iv,ic,exo_nbr);
StateSpaceModel.state_innovations_covariance_matrix = M_.Sigma_e;
StateSpaceModel.order_var = dr.order_var;

conditional_decomposition_array = conditional_variance_decomposition(StateSpaceModel,Steps,SubsetOfVariables );

if options_.noprint == 0
  if options_.order == 2
    disp(' ')                
    disp('CONDITIONAL VARIANCE DECOMPOSITION (in percent), based on first order approximation')
  else
    disp(' ')                
    disp('CONDITIONAL VARIANCE DECOMPOSITION (in percent)')
  end
end

vardec_i = zeros(length(SubsetOfVariables),exo_nbr);

for i=1:length(Steps)
    disp(['Period ' int2str(Steps(i)) ':'])
    
    for j=1:exo_nbr
        vardec_i(:,j) = 100*conditional_decomposition_array(:, ...
                                                          i,j);
    end
    if options_.noprint == 0
        headers = M_.exo_names;
        headers(M_.exo_names_orig_ord,:) = headers;
        headers = char(' ',headers);
        lh = size(deblank(M_.endo_names(SubsetOfVariables,:)),2)+2;
        dyntable('',headers,...
                 deblank(M_.endo_names(SubsetOfVariables,:)),...
                 vardec_i,lh,8,2);
    end
end

oo_.conditional_variance_decomposition = conditional_decomposition_array;
