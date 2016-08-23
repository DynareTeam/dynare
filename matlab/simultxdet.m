function [y_,int_width,int_width_ME]=simultxdet(y0,ex,ex_det, iorder,var_list,M_,oo_,options_)
%function [y_,int_width]=simultxdet(y0,ex,ex_det, iorder,var_list,M_,oo_,options_)
%
% Simulates a stochastic model in the presence of deterministic exogenous shocks
%
% INPUTS:
%    y0:        initial values, of length M_.maximum_lag
%    ex:        matrix of stochastic exogenous shocks, starting at period 1
%    ex_det:    matrix of deterministic exogenous shocks, starting at period 1-M_.maximum_lag
%    iorder:    order of approximation
%    var_list:  list of endogenous variables to simulate
%   int_width_ME:distance between upper bound and
%                mean forecast when considering measurement error
% OUTPUTS:
%   yf:          mean forecast
%   int_width:   distance between upper bound and
%                mean forecast
%   int_width_ME:distance between upper bound and
%                mean forecast when considering measurement error
%
% The forecast horizon is equal to size(ex, 1).
% The condition size(ex,1)+M_.maximum_lag=size(ex_det,1) must be verified
%  for consistency.

% Copyright (C) 2008-2016 Dynare Team
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

dr = oo_.dr;
ykmin = M_.maximum_lag;
endo_nbr = M_.endo_nbr;
exo_det_steady_state = oo_.exo_det_steady_state;
nstatic = M_.nstatic;
nspred = M_.nspred;
nc = size(dr.ghx,2);
iter = size(ex,1);
if size(ex_det, 1) ~= iter+ykmin
    error('Size mismatch: number of forecasting periods for stochastic exogenous and deterministic exogenous don''t match')
end
nx = size(dr.ghu,2);
y_ = zeros(size(y0,1),iter+ykmin);
y_(:,1:ykmin) = y0;
k1 = [ykmin:-1:1];
k2 = dr.kstate(find(dr.kstate(:,2) <= ykmin+1),[1 2]);
k2 = k2(:,1)+(ykmin+1-k2(:,2))*endo_nbr;
k3 = M_.lead_lag_incidence(1:ykmin,:)';
k3 = find(k3(:));
k4 = dr.kstate(find(dr.kstate(:,2) < ykmin+1),[1 2]);
k4 = k4(:,1)+(ykmin+1-k4(:,2))*endo_nbr;

nvar = size(var_list,1);
if nvar == 0
    nvar = endo_nbr;
    ivar = [1:nvar];
else
    ivar=zeros(nvar,1);
    for i=1:nvar
        i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
        if isempty(i_tmp)
            disp(var_list(i,:));
            error (['One of the variable specified does not exist']) ;
        else
            ivar(i) = i_tmp;
        end
    end
end

if iorder == 1
    for i = ykmin+1: iter+ykmin
        tempx1 = y_(dr.order_var,k1);
        tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin);
        tempx = tempx2(k2);
        y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx+dr.ghu* ...
            ex(i-ykmin,:)';
        for j=1:min(ykmin+M_.exo_det_length+1-i,M_.exo_det_length)
            y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*(ex_det(i+j-1,:)'-exo_det_steady_state);
        end
        
        k1 = k1+1;
    end
elseif iorder == 2
    for i = ykmin+1: iter+ykmin
        tempx1 = y_(dr.order_var,k1);
        tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin);
        tempx = tempx2(k2);
        tempu = ex(i-ykmin,:)';
        tempuu = kron(tempu,tempu);
        tempxx = kron(tempx,tempx);
        tempxu = kron(tempx,tempu);
        y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+ ...
            dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu* ...
            tempxu;
        for j=1:min(ykmin+M_.exo_det_length+1-i,M_.exo_det_length)
            tempud = ex_det(i+j-1,:)'-exo_det_steady_state;
            tempudud = kron(tempud,tempud);
            tempxud = kron(tempx,tempud);
            tempuud = kron(tempu,tempud);
            y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*tempud + ...
                dr.ghxud{j}*tempxud + dr.ghuud{j}*tempuud + ...
                0.5*dr.ghudud{j,j}*tempudud;
            for k=1:j-1
                tempudk = ex_det(i+k-1,:)'-exo_det_steady_state;
                tempududk = kron(tempudk,tempud);
                y_(dr.order_var,i) = y_(dr.order_var,i) + ...
                    dr.ghudud{k,j}*tempududk;
            end
        end
        k1 = k1+1;
    end
end

[A,B] = kalman_transition_matrix(dr,nstatic+(1:nspred),1:nc,M_.exo_nbr);

inv_order_var = dr.inv_order_var;
ghx1 = dr.ghx(inv_order_var(ivar),:);
ghu1 = dr.ghu(inv_order_var(ivar),:);

sigma_u = B*M_.Sigma_e*B';
sigma_u1 = ghu1*M_.Sigma_e*ghu1';
sigma_y = 0;

var_yf=NaN(iter,nvar); %initialize
for i=1:iter
    sigma_y1 = ghx1*sigma_y*ghx1'+sigma_u1;
    var_yf(i,:) = diag(sigma_y1)';
    if i == iter
        break
    end
    sigma_u = A*sigma_u*A';
    sigma_y = sigma_y+sigma_u;
end

fact = norminv((1-options_.forecasts.conf_sig)/2,0,1);
if nargout==3
    var_yf_ME=var_yf;
    var_yf_ME(:,options_.varobs_id)=var_yf(:,options_.varobs_id)+repmat(diag(M_.H)',horizon,1);
    int_width_ME = zeros(horizon,M_.endo_nbr);
end

int_width = zeros(iter,nvar);
for i=1:nvar
    int_width(:,i) = fact*sqrt(var_yf(:,i));
        if nargout==3
        int_width_ME(:,i) = -fact*sqrt(var_yf_ME(:,i));
        end
end
