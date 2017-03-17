function [A,B] = ghx2transition(mm,iv,ic,aux)
% [A,B] = ghx2transition(mm,iv,ic,aux)
%
% Adapted by M. Ratto (from kalman_transition_matrix.m)
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 
%

% Copyright (C) 2012 Dynare Team
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

global oo_ M_

  [nr1, nc1] = size(mm);
  ghx = mm(:, [1:(nc1-M_.exo_nbr)]);
  ghu = mm(:, [(nc1-M_.exo_nbr+1):end] );
  if nargin == 1
    oo_.dr.ghx = ghx;
    oo_.dr.ghu = ghu;
    endo_nbr = M_.endo_nbr;
    nstatic = M_.nstatic;
    nspred = M_.nspred;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:nspred) endo_nbr+(1:size(oo_.dr.ghx,2)-nspred) ]';
    aux = oo_.dr.transition_auxiliary_variables;
    k = find(aux(:,2) > nspred);
    aux(:,2) = aux(:,2) + nstatic;
    aux(k,2) = aux(k,2) + M_.nfwrd;
  end
  n_iv = length(iv);
  n_ir1 = size(aux,1);
  nr = n_iv + n_ir1;
  
  A = zeros(nr,nr);
  B = zeros(nr,M_.exo_nbr);
  
  i_n_iv = 1:n_iv;
  A(i_n_iv,ic) = ghx(iv,:);
  if n_ir1 > 0
    A(n_iv+1:end,:) = sparse(aux(:,1),aux(:,2),ones(n_ir1,1),n_ir1,nr);
  end
  
  B(i_n_iv,:) = ghu(iv,:);
