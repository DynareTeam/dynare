function info = endogenous_prior_restrictions(T,R,Model,DynareOptions,DynareResults);
% Check for prior (sign) restrictions on irf's 
%
% INPUTS
%    T          [double]     n*n state space matrix 
%    R          [double]     n*k matrix of shocks
%    Model      [structure]
%    DynareOptions [structure]
%    DynareResults [structure]

% OUTPUTS
%    info     [double]  check if prior restrictions are matched by the
%                       model and related info

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

endo_prior_restrictions= DynareOptions.endogenous_prior_restrictions;
info=0;
if isempty(endo_prior_restrictions.irf),
    % no restriction to be checked
    return
end
if DynareOptions.order>1,
    error('The algorithm for prior (sign) restrictions on irf''s is currently restricted to first-order decision rules')
	return
end
infos=[0 0];
varlist=Model.endo_names(DynareResults.dr.order_var,:);
varlist=varlist(DynareResults.dr.restrict_var_list,:);
T=1;
for j=1:size(endo_prior_restrictions.irf,1),
    T=max(T,endo_prior_restrictions.irf{j,3});
end
for t=1:T,
    RR = T^(t-1)*R;
    for j=1:size(endo_prior_restrictions.irf,1),
	    if endo_prior_restrictions.irf{j,3}~=t,
		   continue,
		end
        iendo=strmatch(endo_prior_restrictions.irf{j,1},varlist,'exact');
        iexo=strmatch(endo_prior_restrictions.irf{j,2},Model.exo_names,'exact');
        if (RR(iendo,iexo)>endo_prior_restrictions.irf{j,4}(1)) && (RR(iendo,iexo)<endo_prior_restrictions.irf{j,4}(2)),
            infos(j,:)=[0, 0];
        else
            if RR(iendo,iexo)<endo_prior_restrictions.irf{j,4}(1),
                delt = (RR(iendo,iexo)-endo_prior_restrictions.irf{j,4}(1))^2;
            else
                delt = (RR(iendo,iexo)-endo_prior_restrictions.irf{j,4}(2))^2;
            end            
            infos(j,:)=[49, delt];
        end
    end
end
if any(infos),
    info=[49,sum(infos(:,2))];
end
return


