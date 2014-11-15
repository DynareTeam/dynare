function T = reduced_rank_cholesky(X)
% Computes the cholesky decomposition of a symetric semidefinite matrix or of a definite positive matrix.

%@info:
%! @deftypefn {Function File} { @var{T} =} reduced_rank_cholesky (@var{X})
%! @anchor{reduced_rank_cholesky}
%! @sp 1
%! Computes the cholesky decomposition of a symetric semidefinite matrix or of a definite positive matrix.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item X
%! n*n matrix of doubles to be factorized (X is supposed to be semidefinite positive).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item T
%! q*n matrix of doubles such that T'*T = X, where q is the number of positive eigenvalues in X.
%! @end table
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! [1] If X is not positive definite, then X has to be a symetric semidefinite matrix.
%! @sp 1
%! [2] The matrix T is upper triangular iff X is positive definite.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{particle/sequential_importance_particle_filter}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2009-2011 Dynare Team
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
    
[T,X_is_not_positive_definite] = chol(X);

if X_is_not_positive_definite
    n = length(X);
    [U,D] = eig(X);
    [tmp,max_elements_indices] = max(abs(U),[],1);
    negloc = (U(max_elements_indices+(0:n:(n-1)*n))<0);
    U(:,negloc) = -U(:,negloc);
    D = diag(D);
    tol = sqrt(eps(max(D))*length(D)*10);
    t = (abs(D) > tol);
    D = D(t);
    if ~(sum(D<0))
        T = diag(sqrt(D))*U(:,t)';
    else
        disp('reduced_rank_cholesky:: Input matrix is not semidefinite positive!')
        T = NaN;
    end
end

%@test:1
%$ n = 10;
%$ m = 100;
%$
%$ X = randn(n,m);
%$ X = X*X';
%$
%$ t = ones(2,1);
%$
%$ try
%$    T = reduced_rank_cholesky(X);
%$ catch
%$    t(1) = 0;
%$    T = all(t);
%$    return;
%$ end
%$
%$
%$ % Check the results.
%$ t(2) = dassert(T,chol(X),1e-16);
%$ T = all(t);
%@eof:1