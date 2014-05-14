function k  = commutation_sparse(n, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% commutation_sparse.m
%
% This file produces the commutation matrix using sparse matrices and was 
% modified from the file commutation.m by Thomas P Minka
% 
%Included in the package nonlinear_MA
%
%THIS VERSION: 1.0.8 Jun. 3, 2013
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% commutation(n, m) or commutation([n m])
% returns Magnus and Neudecker's commutation matrix of dimensions n by m

% Author: Thomas P Minka (tpminka@media.mit.edu)

if nargin < 2
  m = n(2);
  n = n(1);
end

  i = 1:(n*m);
  a = reshape(i, n, m);
  j = vec(a');
  k=sparse(i',j(i),ones(n*m,1));

