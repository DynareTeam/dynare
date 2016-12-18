function tau = thet2tau(params, estim_params_, M_, oo_, indx, indexo, flagmoments,mf,nlags,useautocorr,iv)

%
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

global options_

if nargin==1,
    indx = [1:M_.param_nbr];
    indexo = [];
end

if nargin<7,
    flagmoments=0;
end
if nargin<10 || isempty(useautocorr),
    useautocorr=0;
end
if nargin<11 || isempty(iv),
    iv=[1:M_.endo_nbr];
end

if length(params)>length(indx),
    M_ = set_all_parameters(params,estim_params_,M_);
else
    M_.params(indx) = params;
end
% if ~isempty(indexo)
%     M_.Sigma_e(indexo,indexo) = diag(params(1:length(indexo)).^2);
% end
[A,B,tele,tubbies,M_,options_,oo_] = dynare_resolve(M_,options_,oo_);
if flagmoments==0,
    ys=oo_.dr.ys(oo_.dr.order_var);
    tau = [ys(iv); vec(A(iv,iv)); dyn_vech(B(iv,:)*M_.Sigma_e*B(iv,:)')];
elseif flagmoments==-1
    [I,J]=find(M_.lead_lag_incidence');
    yy0=oo_.dr.ys(I);
    [residual, g1] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1);
    tau=[oo_.dr.ys(oo_.dr.order_var); g1(:)];

else
    GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
    k = find(abs(GAM) < 1e-12);
    GAM(k) = 0;
    if useautocorr,
        sy = sqrt(diag(GAM));
        sy = sy*sy';
        sy0 = sy-diag(diag(sy))+eye(length(sy));
        dum = GAM./sy0;
        tau = dyn_vech(dum(mf,mf));
    else
        tau = dyn_vech(GAM(mf,mf));
    end
    for ii = 1:nlags
        dum = A^(ii)*GAM;
        if useautocorr,
            dum = dum./sy;
        end
        tau = [tau;vec(dum(mf,mf))];
    end
    tau = [ oo_.dr.ys(oo_.dr.order_var(mf)); tau];
end
