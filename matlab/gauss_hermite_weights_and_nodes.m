function [nodes,weights] = gauss_hermite_weights_and_nodes(n)
% Computes the weights and nodes for an Hermite Gaussian quadrature rule.

%@info:
%! @deftypefn {Function File} {@var{nodes}, @var{weights} =} gauss_hermite_weights_and_nodes (@var{n})
%! @anchor{gauss_hermite_weights_and_nodes}
%! @sp 1
%! Computes the weights and nodes for an Hermite Gaussian quadrature rule. designed to approximate integrals
%! on the infinite interval (-\infty,\infty) of an unweighted smooth function.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item n
%! Positive integer scalar, number of nodes (order of approximation).
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item nodes
%! n*1 vector of doubles, the nodes (roots of an order n Hermite polynomial)
%! @item weights
%! n*1 vector of doubles, the associated weights.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2012 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

b = sqrt([1:n-1]/2);
JacobiMatrix = diag(b,1)+diag(b,-1);
[JacobiEigenVectors,JacobiEigenValues] = eig(JacobiMatrix);
[nodes,idx] = sort(diag(JacobiEigenValues));
JacobiEigenVector = JacobiEigenVectors(1,:);
JacobiEigenVector = transpose(JacobiEigenVector(idx));
weights = JacobiEigenVector.^2;
nodes = sqrt(2)*nodes;

%@test:1
%$ n = 5;
%$ [nodes,weights] = gauss_hermite_weights_and_nodes(n);
%$
%$ sum_of_weights = sum(weights);
%$
%$ % Expected nodes (taken from Judd (1998, table 7.4).
%$ enodes = [-2.020182870; -0.9585724646; 0; 0.9585724646;   2.020182870]; 
%$
%$ % Check the results.
%$ t(1) = dassert(1.0,sum_of_weights,1e-12);
%$ t(2) = dassert(enodes,nodes/sqrt(2),1e-8);
%$ T = all(t);
%@eof:1

%@test:2
%$ n = 9;
%$ [nodes,weights] = gauss_hermite_weights_and_nodes(n);
%$
%$ sum_of_weights = sum(weights);
%$ expectation = sum(weights.*nodes); 
%$ variance = sum(weights.*(nodes.^2));
%$
%$ % Check the results.
%$ t(1) = dassert(1.0,sum_of_weights,1e-12);
%$ t(2) = dassert(1.0,variance,1e-12);
%$ t(3) = dassert(0.0,expectation,1e-12);
%$ T = all(t);
%@eof:2

%@test:3
%$ n = 9;
%$ [nodes,weights] = gauss_hermite_weights_and_nodes(n);
%$
%$ NODES = cartesian_product_of_sets(nodes,nodes);
%$ WEIGHTS = cartesian_product_of_sets(weights,weights);
%$ WEIGHTS = prod(WEIGHTS,2);
%$
%$ sum_of_weights = sum(WEIGHTS);
%$ expectation = transpose(WEIGHTS)*NODES;
%$ variance = transpose(WEIGHTS)*NODES.^2;
%$
%$ % Check the results.
%$ t(1) = dassert(1.0,sum_of_weights,1e-12);
%$ t(2) = dassert(ones(1,2),variance,1e-12);
%$ t(3) = dassert(zeros(1,2),expectation,1e-12);
%$ T = all(t);
%@eof:3

%@test:4
%$ n = 9; sigma = .1;
%$ [nodes,weights] = gauss_hermite_weights_and_nodes(n);
%$
%$ sum_of_weights = sum(weights);
%$ expectation = sum(weights.*nodes*.1); 
%$ variance = sum(weights.*((nodes*.1).^2));
%$
%$ % Check the results.
%$ t(1) = dassert(1.0,sum_of_weights,1e-12);
%$ t(2) = dassert(.01,variance,1e-12);
%$ t(3) = dassert(0.0,expectation,1e-12);
%$ T = all(t);
%@eof:4
