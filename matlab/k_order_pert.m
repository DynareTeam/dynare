function [dr,info] = k_order_pert(dr,M,options)
% Compute decision rules using the k-order DLL from Dynare++

% Copyright (C) 2009-2013 Dynare Team
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

info = 0;

M.var_order_endo_names = M.endo_names(dr.order_var,:);

order = options.order;
endo_nbr = M.endo_nbr;
exo_nbr = M.exo_nbr;
nspred = M.nspred;

switch(order)
  case 1
    [err, g_1] = k_order_perturbation(dr,M,options);
    if err
      info(1)=9;
      return;
    end
    dr.g_1 = g_1;
  case 2
    [err, g_0, g_1, g_2] = k_order_perturbation(dr,M,options);
    if err
      info(1)=9;
      return;
    end
    dr.g_0 = g_0;
    dr.g_1 = g_1;
    dr.g_2 = g_2;
  case 3
    if options.pruning 
        [err, g_0, g_1, g_2, g_3, derivs] = k_order_perturbation(dr, ...
                                                          M,options);
        if err
          info(1)=9;
          return;
        end
        dr.ghx = derivs.gy;
        dr.ghu = derivs.gu;
        dr.ghxx = unfold2(derivs.gyy,nspred);
        dr.ghxu = derivs.gyu;
        dr.ghuu = unfold2(derivs.guu,exo_nbr);
        dr.ghs2 = derivs.gss;
        dr.ghxxx = unfold3(derivs.gyyy,nspred);
        dr.ghxxu = unfold21(derivs.gyyu,nspred,exo_nbr);
        dr.ghxuu = unfold12(derivs.gyuu,nspred,exo_nbr);
        dr.ghuuu = unfold3(derivs.guuu,exo_nbr);
        dr.ghxss = derivs.gyss;
        dr.ghuss = derivs.guss;
    else
        [err, g_0, g_1, g_2, g_3] = k_order_perturbation(dr, ...
                                                         M,options);
        if err
          info(1)=9;
          return;
        end
    end
    dr.g_0 = g_0;
    dr.g_1 = g_1;
    dr.g_2 = g_2;
    dr.g_3 = g_3;
  otherwise
    error('order > 3 isn''t implemented')
end

if options.pruning 
    return
end

nspred = M.nspred;

dr.ghx = dr.g_1(:,1:nspred);
dr.ghu = dr.g_1(:,nspred+1:end);

if options.loglinear == 1
    k = find(dr.kstate(:,2) <= M.maximum_endo_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;
    
    dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
             repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
    dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
end

if order > 1
    dr.ghs2 = 2*g_0;
    s0 = 0;
    s1 = 0;
    ghxx=zeros(endo_nbr, nspred^2);
    ghxu=zeros(endo_nbr, nspred*exo_nbr);
    ghuu=zeros(endo_nbr, exo_nbr^2);
    for i=1:size(g_2,2)
        if s0 < nspred && s1 < nspred
            ghxx(:,s0*nspred+s1+1) = 2*g_2(:,i);
            if s1 > s0
                ghxx(:,s1*nspred+s0+1) = 2*g_2(:,i);
            end
        elseif s0 < nspred && s1 < nspred+exo_nbr 
            ghxu(:,(s0*exo_nbr+s1-nspred+1)) = 2*g_2(:,i);
        elseif s0 < nspred+exo_nbr && s1 < nspred+exo_nbr
            ghuu(:,(s0-nspred)*exo_nbr+s1-nspred +1) = 2*g_2(:,i);
            if s1 > s0
                ghuu(:,(s1-nspred)*exo_nbr+s0-nspred+1) = 2*g_2(:,i);
            end
        else
            error('dr1:k_order_perturbation:g_2','Unaccounted columns in g_2');
        end
        s1 = s1+1;
        if s1 == nspred+exo_nbr
            s0 = s0+1;
            s1 = s0; 
        end
    end % for loop
    dr.ghxx = ghxx;
    dr.ghxu = ghxu;
    dr.ghuu = ghuu;
end

function y = unfold2(x,n)
y=zeros(size(x,1),n*n);
m = 1;
for i=1:n
    for j=i:n
        y(:,(i-1)*n+j)=x(:,m);
        if j ~= i
            y(:,(j-1)*n+i)=x(:,m);
        end
        m = m+1;
    end
end

function y = unfold3(x,n)
y = zeros(size(x,1),n*n*n);
m = 1;
for i=1:n
    for j=i:n
        for k=j:n
            xx = x(:,m);
            y(:,(i-1)*n*n+(j-1)*n+k) = xx;
            y(:,(i-1)*n*n+(k-1)*n+j) = xx;
            y(:,(j-1)*n*n+(k-1)*n+i) = xx;
            y(:,(j-1)*n*n+(i-1)*n+k) = xx;
            y(:,(k-1)*n*n+(i-1)*n+j) = xx;
            y(:,(k-1)*n*n+(j-1)*n+i) = xx;
            m = m + 1;
        end
    end
end

function y = unfold21(x,n1,n2)
y = zeros(size(x,1),n1*n1*n2);
m = 1;
for i=1:n1
    for j=i:n1
        for k=1:n2
            xx = x(:,m);
            y(:,(i-1)*n1*n2+(j-1)*n2+k) = xx;
            if j ~= i
                y(:,(j-1)*n1*n2+(i-1)*n2+k) = xx;
            end
            m = m + 1;
        end
    end
end 

function y = unfold12(x,n1,n2)
y = zeros(size(x,1),n1*n2*n2);
m = 1;
for i=1:n1
    for j=1:n2
        for k=j:n2
            xx = x(:,m);
            y(:,(i-1)*n2*n2+(j-1)*n2+k) = xx;
            if k ~= j
                y(:,(i-1)*n2*n2+(k-1)*n2+j) = xx;
            end
            m = m + 1;
        end
    end
end


            
            
