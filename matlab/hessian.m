function hessian_mat = hessian(func,x,gstep,varargin) % --*-- Unitary tests --*--

% Computes second order partial derivatives
%
% INPUTS
%    func        [string]   name of the function
%    x           [double]   vector, the Hessian of "func" is evaluated at x.
%    gstep       [double]   scalar, size of epsilon.
%    varargin    [void]     list of additional arguments for "func".
%
% OUTPUTS
%    hessian_mat [double]   Hessian matrix
%
% ALGORITHM
%    Uses Abramowitz and Stegun (1965) formulas 25.3.23 
% \[
%     \frac{\partial^2 f_{0,0}}{\partial {x^2}} = \frac{1}{h^2}\left( f_{1,0} - 2f_{0,0} + f_{ - 1,0} \right)
% \]
% and 25.3.27 p. 884
% 
% \[
%     \frac{\partial ^2f_{0,0}}{\partial x\partial y} = \frac{-1}{2h^2}\left(f_{1,0} + f_{-1,0} + f_{0,1} + f_{0,-1} - 2f_{0,0} - f_{1,1} - f_{-1,-1} \right)
% \]
%
% SPECIAL REQUIREMENTS
%    none
%  

% Copyright (C) 2001-2014 Dynare Team
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

if ~isa(func, 'function_handle') 
    func = str2func(func);
end
n=size(x,1);
h1=max(abs(x),sqrt(gstep(1))*ones(n,1))*eps^(1/6)*gstep(2);
h_1=h1;
xh1=x+h1;
h1=xh1-x;
xh1=x-h_1;
h_1=x-xh1;
xh1=x;
f0=feval(func,x,varargin{:});
f1=zeros(size(f0,1),n);
f_1=f1;
for i=1:n
    %do step up
    xh1(i)=x(i)+h1(i);
    f1(:,i)=feval(func,xh1,varargin{:});
    %do step up
    xh1(i)=x(i)-h_1(i);
    f_1(:,i)=feval(func,xh1,varargin{:});
    xh1(i)=x(i);%reset parameter
end
xh_1=xh1;
hessian_mat = zeros(size(f0,1),n*n);
temp=f1+f_1-f0*ones(1,n); %term f_(1,0)+f_(-1,0)-f_(0,0) used later
for i=1:n    
    if i > 1  %fill symmetric part of Hessian based on previously computed results      
        k=[i:n:n*(i-1)];
        hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
    end     
    hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i)); %formula 25.3.23
    for j=i+1:n        
        %step in up direction
        xh1(i)=x(i)+h1(i);
        xh1(j)=x(j)+h_1(j);
        %step in down direction
        xh_1(i)=x(i)-h1(i);
        xh_1(j)=x(j)-h_1(j);
        hessian_mat(:,(i-1)*n+j)=-(-feval(func,xh1,varargin{:})-feval(func,xh_1,varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j)); %formula 25.3.27
        %reset grid points
        xh1(i)=x(i);
        xh1(j)=x(j);
        xh_1(i)=x(i);
        xh_1(j)=x(j);
    end    
end


%@test:1
%$ % Create a function.
%$ fid = fopen('exfun.m','w+');
%$ fprintf(fid,'function [f,g,H] = exfun(xvar)\\n');
%$ fprintf(fid,'x = xvar(1);\\n');
%$ fprintf(fid,'y = xvar(2);\\n');
%$ fprintf(fid,'f = x^2* log(y);\\n');
%$ fprintf(fid,'if nargout>1\\n');
%$ fprintf(fid,'    g = zeros(2,1);\\n');
%$ fprintf(fid,'    g(1) = 2*x*log(y);\\n');
%$ fprintf(fid,'    g(2) = x*x/y;\\n');
%$ fprintf(fid,'end\\n');
%$ fprintf(fid,'if nargout>2\\n');
%$ fprintf(fid,'    H = zeros(2,2);\\n');
%$ fprintf(fid,'    H(1,1) = 2*log(y);\\n');
%$ fprintf(fid,'    H(1,2) = 2*x/y;\\n');
%$ fprintf(fid,'    H(2,1) = H(1,2);\\n');
%$ fprintf(fid,'    H(2,2) = -x*x/(y*y);\\n');
%$ fprintf(fid,'    H = H(:);\\n');
%$ fprintf(fid,'end\\n');
%$ fclose(fid);
%$
%$ rehash;
%$
%$ t = zeros(5,1);
%$
%$ % Evaluate the Hessian at (1,e)
%$ try
%$    H = hessian('exfun',[1; exp(1)],[1e-2; 1]);
%$    t(1) = 1;
%$ catch
%$    t(1) = 0;
%$ end
%$
%$ % Compute the true Hessian matrix
%$ [f, g, Htrue] = exfun([1 exp(1)]);
%$
%$ % Delete exfun routine from disk.
%$ delete('exfun.m');
%$
%$ % Compare the values in H and Htrue
%$ if t(1)
%$    t(2) = dassert(abs(H(1)-Htrue(1))<1e-6,true);
%$    t(3) = dassert(abs(H(2)-Htrue(2))<1e-6,true);
%$    t(4) = dassert(abs(H(3)-Htrue(3))<1e-6,true);
%$    t(5) = dassert(abs(H(4)-Htrue(4))<1e-6,true);
%$ end
%$ T = all(t);
%@eof:1
