function [X, info] = cycle_reduction(A0, A1, A2, cvg_tol, ch) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {[@var{X}, @var{info}] =} cycle_reduction (@var{A0},@var{A1},@var{A2},@var{cvg_tol},@var{ch})
%! @anchor{cycle_reduction}
%! @sp 1
%! Solves the quadratic matrix equation A2*X^2 + A1*X + A0 = 0.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A0
%! Square matrix of doubles, n*n.
%! @item A1
%! Square matrix of doubles, n*n.
%! @item A2
%! Square matrix of doubles, n*n.
%! @item cvg_tol
%! Scalar double, tolerance parameter.
%! @item ch
%! Any matlab object, if not empty the solution is checked.
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
%! @sp 2
%! @strong{References:}
%! @sp 1
%! D.A. Bini, G. Latouche, B. Meini (2002), "Solving matrix polynomial equations arising in queueing problems", Linear Algebra and its Applications 340, pp. 222-244
%! @sp 1
%! D.A. Bini, B. Meini (1996), "On the solution of a nonlinear matrix equation arising in queueing problems", SIAM J. Matrix Anal. Appl. 17, pp. 906-926.
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2013 Dynare Team
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

max_it = 300;
it = 0;
info = 0;
X = [];
crit = Inf;
A0_0 = A0;
Ahat1 = A1;
if (nargin == 5 && ~isempty(ch) )
    A1_0 = A1;
    A2_0 = A2;
end
n = length(A0);
id0 = 1:n;
id2 = id0+n;

cont = 1;
while cont 
    tmp = ([A0; A2]/A1)*[A0 A2];
    A1 = A1 - tmp(id0,id2) - tmp(id2,id0);
    A0 = -tmp(id0,id0);
    A2 = -tmp(id2,id2);
    Ahat1 = Ahat1 -tmp(id2,id0);
    crit = norm(A0,1);
    if crit < cvg_tol
        % keep iterating until condition on A2 is met
        if norm(A2,1) < cvg_tol
            cont = 0;
        end
    elseif isnan(crit) || it == max_it
        if crit < cvg_tol
            info(1) = 4;
            info(2) = log(norm(A2,1));
        else
            info(1) = 3;
            info(2) = log(norm(A1,1));
        end
        return
    end        
    it = it + 1;
end

X = -Ahat1\A0_0;

if (nargin == 5 && ~isempty(ch) )
    %check the solution
    res = A0_0 + A1_0 * X + A2_0 * X * X;
    if (sum(sum(abs(res))) > cvg_tol)
        disp(['the norm residual of the residu ' num2str(res) ' compare to the tolerance criterion ' num2str(cvg_tol)]);
    end
end

%@test:1
%$
%$ t = zeros(3,1);
%$
%$ % Set the dimension of the problem to be solved.
%$ n = 500;
%$
%$ % Set the equation to be solved
%$ A = eye(n);
%$ B = diag(30*ones(n,1)); B(1,1) = 20; B(end,end) = 20; B = B - diag(10*ones(n-1,1),-1); B = B - diag(10*ones(n-1,1),1);
%$ C = diag(15*ones(n,1)); C = C - diag(5*ones(n-1,1),-1); C = C - diag(5*ones(n-1,1),1);
%$
%$ % Solve the equation with the cycle reduction algorithm
%$ try
%$     t=cputime; X1 = cycle_reduction(C,B,A,1e-7); elapsedtime = cputime-t;
%$     disp(['cputime for cycle reduction algorithm is: ' num2str(elapsedtime) ' (n=' int2str(n) ').'])
%$     t(1) = 1;
%$ catch
%$     % nothing to do here.
%$ end
%$
%$ % Solve the equation with the logarithmic reduction algorithm
%$ try
%$     t=cputime; X2 = logarithmic_reduction(A,B,C,1e-16,100); elapsedtime = cputime-t;
%$     disp(['cputime for logarithmic reduction algorithm is: ' num2str(elapsedtime) ' (n=' int2str(n) ').'])
%$     t(2) = 1;
%$ catch
%$     % nothing to do here.
%$ end
%$
%$ % Check the results.
%$ if t(1) && t(2)
%$     t(3) = dassert(X1,X2,1e-12);
%$ end
%$
%$ T = all(t);
%@eof:1