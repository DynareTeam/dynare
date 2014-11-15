function [nodes,weights] = gauss_legendre_weights_and_nodes(n,a,b)
% Computes the weights and nodes for a Legendre Gaussian quadrature rule.

%@info:
%! @deftypefn {Function File} {@var{nodes}, @var{weights} =} gauss_hermite_weights_and_nodes (@var{n})
%! @anchor{gauss_legendre_weights_and_nodes}
%! @sp 1
%! Computes the weights and nodes for a Legendre Gaussian quadrature rule. designed to approximate integrals
%! on the finite interval (-1,1) of an unweighted smooth function.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item n
%! Positive integer scalar, number of nodes (order of approximation).
%! @item a
%! Double scalar, lower bound.
%! @item b
%! Double scalar, upper bound.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item nodes
%! n*1 vector of doubles, the nodes (roots of an order n Legendre polynomial)
%! @item weights
%! n*1 vector of doubles, the associated weights.
%! @end table
%! @sp 2
%! @strong{Remarks:}
%! Only the first input argument (the number of nodes) is mandatory. The second and third input arguments
%! are used if a change of variables is necessary (ie if we need nodes over the interval [a,b] instead of
%! of the default interval [-1,1]).
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

bb = sqrt(1./(4-(1./transpose(1:n-1)).^2));
aa = zeros(n,1);
JacobiMatrix = diag(bb,1)+diag(aa)+diag(bb,-1);
[JacobiEigenVectors,JacobiEigenValues] = eig(JacobiMatrix);
[nodes,idx] = sort(diag(JacobiEigenValues));
JacobiEigenVector = JacobiEigenVectors(1,:);
JacobiEigenVector = transpose(JacobiEigenVector(idx));
weights = 2*JacobiEigenVector.^2;

if nargin==3
    weights = .5*(b-a)*weights;
    nodes = .5*(nodes+1)*(b-a)+a;
end

%@test:1
%$ [n2,w2] = gauss_legendre_weights_and_nodes(2);
%$ [n3,w3] = gauss_legendre_weights_and_nodes(3);
%$ [n4,w4] = gauss_legendre_weights_and_nodes(4);
%$ [n5,w5] = gauss_legendre_weights_and_nodes(5);
%$ [n7,w7] = gauss_legendre_weights_and_nodes(7);
%$
%$
%$ % Expected nodes (taken from Judd (1998, table 7.2)).
%$ e2 = .5773502691; e2 = [-e2; e2];
%$ e3 = .7745966692; e3 = [-e3; 0 ; e3];
%$ e4 = [.8611363115; .3399810435]; e4 = [-e4; flipud(e4)];
%$ e5 = [.9061798459; .5384693101]; e5 = [-e5; 0; flipud(e5)];
%$ e7 = [.9491079123; .7415311855; .4058451513]; e7 = [-e7; 0; flipud(e7)];
%$
%$ % Expected weights (taken from Judd (1998, table 7.2) and http://en.wikipedia.org/wiki/Gaussian_quadrature).
%$ f2 = [1; 1];
%$ f3 = [5; 8; 5]/9;
%$ f4 = [18-sqrt(30); 18+sqrt(30)]; f4 = [f4; flipud(f4)]/36;
%$ f5 = [322-13*sqrt(70); 322+13*sqrt(70)]/900; f5 = [f5; 128/225; flipud(f5)];
%$ f7 = [.1294849661; .2797053914; .3818300505]; f7 = [f7; .4179591836; flipud(f7)];
%$
%$ % Check the results.
%$ t(1) = dassert(e2,n2,1e-9);
%$ t(2) = dassert(e3,n3,1e-9);
%$ t(3) = dassert(e4,n4,1e-9);
%$ t(4) = dassert(e5,n5,1e-9);
%$ t(5) = dassert(e7,n7,1e-9);
%$ t(6) = dassert(w2,f2,1e-9);
%$ t(7) = dassert(w3,f3,1e-9);
%$ t(8) = dassert(w4,f4,1e-9);
%$ t(9) = dassert(w5,f5,1e-9);
%$ t(10) = dassert(w7,f7,1e-9);
%$ T = all(t);
%@eof:1

%@test:2
%$ nmax = 50;
%$
%$ t = zeros(nmax,1);
%$
%$ for i=1:nmax
%$     [n,w] = gauss_legendre_weights_and_nodes(i);
%$     t(i) = dassert(sum(w),2,1e-12);
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ [n,w] = gauss_legendre_weights_and_nodes(9,pi,2*pi);
%$ % Check that the 
%$ t(1) = all(n>pi);
%$ t(2) = all(n<2*pi);
%$ t(3) = dassert(sum(w),pi,1e-12);
%$ T = all(t);
%@eof:3