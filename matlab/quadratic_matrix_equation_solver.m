function [X,info] = quadratic_matrix_equation_solver(A,B,C,tol,maxit,line_search_flag,X)

%@info:
%! @deftypefn {Function File} {[@var{X1}, @var{info}] =} quadratic_matrix_equation_solver (@var{A},@var{B},@var{C},@var{tol},@var{maxit},@var{line_search_flag},@var{X0})
%! @anchor{logarithmic_reduction}
%! @sp 1
%! Solves the quadratic matrix equation AX^2 + BX + C = 0 with a Newton algorithm.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Square matrix of doubles, n*n.
%! @item B
%! Square matrix of doubles, n*n.
%! @item C
%! Square matrix of doubles, n*n.
%! @item tol
%! Scalar double, tolerance parameter.
%! @item maxit
%! Scalar integer, maximum number of iterations.
%! @item line_search_flag
%! Scalar integer, if nonzero an exact line search algorithm is used.
%! @item X
%! Square matrix of doubles, n*n, initial condition.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item X
%! Square matrix of doubles, n*n, solution of the matrix equation.
%! @item info
%! Scalar integer, if nonzero the algorithm failed in finding the solution of the matrix equation.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{fastgensylv}
%! @sp 2
%! @strong{References:}
%! @sp 1
%! N.J. Higham and H.-M. Kim (2001), "Solving a quadratic matrix equation by Newton's method with exact line searches.", in SIAM J. Matrix Anal. Appl., Vol. 23, No. 3, pp. 303-316.
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

provide_initial_condition_to_fastgensylv = 0;

info = 0;

F = eval_quadratic_matrix_equation(A,B,C,X);

if max(max(abs(F)))<tol
    return
end

kk = 0.0;
cc = 1+tol;

step_length = 1.0;

while kk<maxit && cc>tol
    if provide_initial_condition_to_fastgensylv && exist('H','var')
        H = fastgensylv(A*X+B,A,X,F,tol,maxit,H);
    else
        try
            H = fastgensylv(A*X+B,A,X,F,tol,maxit);
        catch
            X = zeros(length(X));
            H = fastgensylv(A*X+B,A,X,F,tol,maxit);
        end
    end
    if line_search_flag
        step_length = line_search(A,H,F);
    end
    X = X + step_length*H;
    F = eval_quadratic_matrix_equation(A,B,C,X);
    cc = max(max(abs(F)));
    kk = kk +1;
end

if cc>tol
    X = NaN(size(X));
    info = 1;
end


function f = eval_quadratic_matrix_equation(A,B,C,X)
    f = C + (B + A*X)*X;

function [p0,p1] = merit_polynomial(A,H,F)
    AHH = A*H*H;
    gamma = norm(AHH,'fro')^2;
    alpha = norm(F,'fro')^2;
    beta  = trace(F*AHH*AHH*F);
    p0 = [gamma, -beta, alpha+beta, -2*alpha, alpha];
    p1 = [4*gamma, -3*beta, 2*(alpha+beta), -2*alpha];

function t = line_search(A,H,F)
    [p0,p1] = merit_polynomial(A,H,F);
    if any(isnan(p0)) || any(isinf(p0))
        t = 1.0;
        return
    end
    r = roots(p1);
    s = [Inf(3,1),r];
    for i = 1:3
        if isreal(r(i))
            s(i,1) = p0(1)*r(i)^4 + p0(2)*r(i)^3 + p0(3)*r(i)^2 + p0(4)*r(i) + p0(5);
        end
    end
    s = sortrows(s,1);
    t = s(1,2);
    if t<=1e-12 || t>=2
        t = 1;
    end

%@test:1
%$ addpath ../matlab
%$
%$ % Set the dimension of the problem to be solved
%$ n = 200;
%$ % Set the equation to be solved
%$ A = eye(n);
%$ B = diag(30*ones(n,1)); B(1,1) = 20; B(end,end) = 20; B = B - diag(10*ones(n-1,1),-1); B = B - diag(10*ones(n-1,1),1);
%$ C = diag(15*ones(n,1)); C = C - diag(5*ones(n-1,1),-1); C = C - diag(5*ones(n-1,1),1);
%$
%$ % Solve the equation with the cycle reduction algorithm
%$ tic, X1 = cycle_reduction(C,B,A,1e-7); toc
%$
%$ % Solve the equation with the logarithmic reduction algorithm
%$ tic, X2 = quadratic_matrix_equation_solver(A,B,C,1e-16,100,1,zeros(n)); toc
%$
%$ % Check the results.
%$ t(1) = dassert(X1,X2,1e-12);
%$
%$ T = all(t);
%@eof:1