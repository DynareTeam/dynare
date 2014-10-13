function [x,u] = lyapunov_symm(a,b,third_argument,lyapunov_complex_threshold,method, R)
% Solves the Lyapunov equation x-a*x*a' = b, for b and x symmetric matrices.
% If a has some unit roots, the function computes only the solution of the stable subsystem.
%  
% INPUTS:
%   a                           [double]    n*n matrix.
%   b                           [double]    n*n matrix.
%   third_argument              [double]    scalar, if method <= 2 :
%                                                      qz_criterium = third_argument unit root threshold for eigenvalues in a,
%                                                    elseif method = 3 :
%                                                      tol =third_argument the convergence criteria for fixed_point algorithm.
%   lyapunov_complex_threshold  [double]    scalar, complex block threshold for the upper triangular matrix T.
%   method                      [integer]   Scalar, if method=0 [default] then U, T, n and k are not persistent.  
%                                                      method=1 then U, T, n and k are declared as persistent 
%                                                               variables and the Schur decomposition is triggered.    
%                                                      method=2 then U, T, n and k are declared as persistent 
%                                                               variables and the Schur decomposition is not performed.
%                                                      method=3 fixed point method
% OUTPUTS
%   x:      [double]    m*m solution matrix of the lyapunov equation, where m is the dimension of the stable subsystem.
%   u:      [double]    Schur vectors associated with unit roots  
%
% ALGORITHM
%   Uses reordered Schur decomposition (Bartels-Stewart algorithm)
%   [method<3] or a fixed point algorithm (method==4)
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2006-2012 Dynare Team
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
if nargin<5
    method = 0;
end

if method == 3
    persistent X method1;
    if ~isempty(method1)
        method = method1;
    end;
    tol = third_argument;
%     fprintf('lyapunov_symm:: [method=%d] \n',method);
    if method == 3
        %tol = 1e-10;
        it_fp = 0;
        evol = 100;
        if isempty(X) || length(X)~=length(b)
            X = b;
            max_it_fp = 2000;
        else
            max_it_fp = 300;
        end;
        at = a';
        %fixed point iterations
        while evol > tol && it_fp < max_it_fp;
            X_old = X;
            X = a * X * at + b;
            evol = max(sum(abs(X - X_old))); %norm_1
            %evol = max(sum(abs(X - X_old)')); %norm_inf
            it_fp = it_fp + 1;
        end;
%         fprintf('lyapunov_symm:: lyapunov fixed_point iterations=%d norm=%g\n',it_fp,evol);
        if it_fp >= max_it_fp
            disp(['lyapunov_symm:: convergence not achieved in solution of Lyapunov equation after ' int2str(it_fp) ' iterations, switching method from 3 to 0']);
            method1 = 0;
            method = 0;
        else
            method1 = 3;
            x = X;
            return;
        end;
    end;
elseif method == 4
    % works only with Matlab System Control toolbox or octave the control package,
    if isoctave
        if ~user_has_octave_forge_package('control')
            error('lyapunov=square_root_solver not available; you must install the control package from Octave Forge')
        end
    else
        if ~user_has_matlab_license('control_toolbox')
            error('lyapunov=square_root_solver not available; you must install the control system toolbox')
        end
    end
    chol_b = R*chol(b,'lower');
    Rx = dlyapchol(a,chol_b);
    x = Rx' * Rx;
    return;
end;

qz_criterium = third_argument;
if method
    persistent U T k n
else
    %    if exist('U','var')
    %        clear('U','T','k','n')
    %    end
end

u = [];

if size(a,1) == 1
    x=b/(1-a*a);
    return
end

if method<2
    [U,T] = schur(a);
    e1 = abs(ordeig(T)) > 2-qz_criterium;
    k = sum(e1);       % Number of unit roots. 
    n = length(e1)-k;  % Number of stationary variables.
    if k > 0
        % Selects stable roots
        [U,T] = ordschur(U,T,e1);
        T = T(k+1:end,k+1:end);
    end
end

B = U(:,k+1:end)'*b*U(:,k+1:end);

x = zeros(n,n);
i = n;

while i >= 2
    if abs(T(i,i-1))<lyapunov_complex_threshold
        if i == n
            c = zeros(n,1);
        else
            c = T(1:i,:)*(x(:,i+1:end)*T(i,i+1:end)') + ...
                T(i,i)*T(1:i,i+1:end)*x(i+1:end,i);
        end
        q = eye(i)-T(1:i,1:i)*T(i,i);
        x(1:i,i) = q\(B(1:i,i)+c);
        x(i,1:i-1) = x(1:i-1,i)';
        i = i - 1;
    else
        if i == n
            c = zeros(n,1);
            c1 = zeros(n,1);
        else
            c = T(1:i,:)*(x(:,i+1:end)*T(i,i+1:end)') + ...
                T(i,i)*T(1:i,i+1:end)*x(i+1:end,i) + ...
                T(i,i-1)*T(1:i,i+1:end)*x(i+1:end,i-1);
            c1 = T(1:i,:)*(x(:,i+1:end)*T(i-1,i+1:end)') + ...
                 T(i-1,i-1)*T(1:i,i+1:end)*x(i+1:end,i-1) + ...
                 T(i-1,i)*T(1:i,i+1:end)*x(i+1:end,i);
        end
        q = [  eye(i)-T(1:i,1:i)*T(i,i) ,  -T(1:i,1:i)*T(i,i-1) ; ...
               -T(1:i,1:i)*T(i-1,i)     ,   eye(i)-T(1:i,1:i)*T(i-1,i-1) ];
        z =  q\[ B(1:i,i)+c ; B(1:i,i-1) + c1 ];
        x(1:i,i) = z(1:i);
        x(1:i,i-1) = z(i+1:end);
        x(i,1:i-1) = x(1:i-1,i)';
        x(i-1,1:i-2) = x(1:i-2,i-1)';
        i = i - 2;
    end
end
if i == 1
    c = T(1,:)*(x(:,2:end)*T(1,2:end)') + T(1,1)*T(1,2:end)*x(2:end,1);
    x(1,1) = (B(1,1)+c)/(1-T(1,1)*T(1,1));
end
x = U(:,k+1:end)*x*U(:,k+1:end)';
u = U(:,1:k);
