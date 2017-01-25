function vx1 = get_variance_of_endogenous_variables(dr,i_var)

% function vx1 = get_variance_of_endogenous_variables(dr,i_var)
% Gets the variance of a variables subset
%
% INPUTS
%    dr:        structure of decisions rules for stochastic simulations
%    i_var:     indices of a variables list
%        
% OUTPUTS
%    vx1:       variance-covariance matrix
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2012 Dynare Team
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

global M_ options_

endo_nbr = M_.endo_nbr;

Sigma_e = M_.Sigma_e;

nstatic = M_.nstatic;
nspred = M_.nspred;
ghx = dr.ghx(i_var,:);
ghu = dr.ghu(i_var,:);
nc = size(ghx,2);
n = length(i_var);

[A,B] = kalman_transition_matrix(dr,nstatic+(1:nspred),1:nc,M_.exo_nbr);

[vx,u] = lyapunov_symm(A,B*Sigma_e*B',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold, [], options_.debug);

if size(u,2) > 0
    i_stat = find(any(abs(ghx*u) < options_.Schur_vec_tol,2)); %only set those variances of objective function for which variance is finite     
    ghx = ghx(i_stat,:);
    ghu = ghu(i_stat,:);
else
    i_stat = (1:n)';
end

vx1 = Inf*ones(n,n);
vx1(i_stat,i_stat) = ghx*vx*ghx'+ghu*Sigma_e*ghu';

