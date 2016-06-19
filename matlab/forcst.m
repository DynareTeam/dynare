function [yf,int_width,int_width_ME]=forcst(dr,y0,horizon,var_list,M_,oo_,options_)
% function [yf,int_width,int_width_ME]=forecst(dr,y0,horizon,var_list,M_,oo_,options_)
%   computes mean forecast for a given value of the parameters
%   computes also confidence band for the forecast    
%
% INPUTS:
%   dr:          structure containing decision rules
%   y0:          initial values
%   horizon:     nbr of periods to forecast
%   var_list:    list of variables (character matrix)
%   M_:          Dynare model structure
%   options_:    Dynare options structure
%   oo_:         Dynare results structure

% OUTPUTS:
%   yf:          mean forecast
%   int_width:   distance between upper bound and
%                mean forecast
%   int_width_ME:distance between upper bound and
%                mean forecast when considering measurement error
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2016 Dynare Team
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

oo_=make_ex_(M_,options_,oo_);
yf = simult_(y0,dr,zeros(horizon,M_.exo_nbr),1);
nstatic = M_.nstatic;
nspred = M_.nspred;
nc = size(dr.ghx,2);
endo_nbr = M_.endo_nbr;
inv_order_var = dr.inv_order_var;
[A,B] = kalman_transition_matrix(dr,nstatic+(1:nspred),1:nc,M_.exo_nbr);

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr,:);
end
nvar = size(var_list,1);
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

ghx1 = dr.ghx(inv_order_var(ivar),:);
ghu1 = dr.ghu(inv_order_var(ivar),:);

%initialize recursion
sigma_u = B*M_.Sigma_e*B'; 
sigma_u1 = ghu1*M_.Sigma_e*ghu1';
sigma_y = 0; %no uncertainty about the states

var_yf=NaN(horizon,nvar); %initialize
for i=1:horizon
    %map uncertainty about states into uncertainty about observables
    sigma_y1 = ghx1*sigma_y*ghx1'+sigma_u1;
    var_yf(i,:) = diag(sigma_y1)';
    if i == horizon
        break
    end
    %update uncertainty about states
    sigma_u = A*sigma_u*A';
    sigma_y = sigma_y+sigma_u;
end
if nargout==3
    var_yf_ME=var_yf;
    [loc_H,loc_varlist]=ismember(options_.varobs',options_.varlist);
    loc_varlist(loc_varlist==0)=[];    
    var_yf_ME(:,loc_varlist)=var_yf(:,loc_varlist)+repmat(diag(M_.H(loc_H,loc_H))',horizon,1);
    int_width_ME = zeros(horizon,nvar);
end

fact = norminv((1-options_.forecasts.conf_sig)/2,0,1);

int_width = zeros(horizon,nvar);
for i=1:nvar
    int_width(:,i) = -fact*sqrt(var_yf(:,i));
    if nargout==3
        int_width_ME(:,i) = -fact*sqrt(var_yf_ME(:,i));
    end
end

yf = yf(ivar,:);
