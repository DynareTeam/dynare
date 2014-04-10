function disp_th_moments(dr,var_list)
% Display theoretical moments of variables

% Copyright (C) 2001-2013 Dynare Team
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

global M_ oo_ options_

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
end
nvar = size(var_list,1);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
    if isempty(i_tmp)
        error (['One of the variable specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

[oo_.gamma_y,stationary_vars] = th_autocovariances(dr,ivar,M_,options_);
m = dr.ys(ivar);
non_stationary_vars = setdiff(1:length(ivar),stationary_vars);
m(non_stationary_vars) = NaN;

i1 = find(abs(diag(oo_.gamma_y{1})) > 1e-12);
s2 = diag(oo_.gamma_y{1});
sd = sqrt(s2);
if options_.order == 2
    m = m+oo_.gamma_y{options_.ar+3};
end

z = [ m sd s2 ];
oo_.mean = m;
oo_.var = oo_.gamma_y{1};

if M_.exo_nbr > 1 && size(stationary_vars, 1) > 0
    oo_.variance_decomposition=100*oo_.gamma_y{options_.ar+2};
    if ~options_.noprint %options_.nomoments == 0
        if options_.order == 2
            title='APROXIMATED THEORETICAL MOMENTS';
        else
            title='THEORETICAL MOMENTS';
        end
        if options_.hp_filter
            title = [title ' (HP filter, lambda = ' num2str(options_.hp_filter) ')'];
        end
        headers=char('VARIABLE','MEAN','STD. DEV.','VARIANCE');
        labels = deblank(M_.endo_names(ivar,:));
        lh = size(labels,2)+2;
        dyntable(title,headers,labels,z,lh,11,4);

        skipline()
            if options_.order == 2
                title='APPROXIMATED VARIANCE DECOMPOSITION (in percent)';            
            else
                title='VARIANCE DECOMPOSITION (in percent)';
            end
            if options_.hp_filter
                title = [title ' (HP filter, lambda = ' ...
                         num2str(options_.hp_filter) ')'];
            end
            headers = M_.exo_names;
            headers(M_.exo_names_orig_ord,:) = headers;
            headers = char(' ',headers);
            lh = size(deblank(M_.endo_names(ivar(stationary_vars),:)),2)+2;
            dyntable(title,headers,deblank(M_.endo_names(ivar(stationary_vars), ...
                                                         :)),100*oo_.gamma_y{options_.ar+2}(stationary_vars,:),lh,8,2);
    end
    
    conditional_variance_steps = options_.conditional_variance_decomposition;
    if length(conditional_variance_steps)
        StateSpaceModel.number_of_state_equations = M_.endo_nbr;
        StateSpaceModel.number_of_state_innovations = M_.exo_nbr;
        StateSpaceModel.sigma_e_is_diagonal = M_.sigma_e_is_diagonal;
        [StateSpaceModel.transition_matrix,StateSpaceModel.impulse_matrix] = kalman_transition_matrix(dr,(1:M_.endo_nbr)',M_.nstatic+(1:M_.nspred)',M_.exo_nbr);
        StateSpaceModel.state_innovations_covariance_matrix = M_.Sigma_e;
        StateSpaceModel.order_var = dr.order_var;
        oo_.conditional_variance_decomposition = conditional_variance_decomposition(StateSpaceModel,conditional_variance_steps,ivar);
        
        if options_.noprint == 0
            display_conditional_variance_decomposition(oo_.conditional_variance_decomposition,conditional_variance_steps,...
                                                         ivar,M_,options_);
        end
    end
end

if length(i1) == 0
    skipline()
    disp('All endogenous are constant or non stationary, not displaying correlations and auto-correlations')
    skipline()
    return
end

if options_.nocorr == 0 && size(stationary_vars, 1) > 0
    corr = oo_.gamma_y{1}(i1,i1)./(sd(i1)*sd(i1)');
    if ~options_.noprint,
        skipline()
        if options_.order == 2
            title='APPROXIMATED MATRIX OF CORRELATIONS';            
        else
            title='MATRIX OF CORRELATIONS';
        end
        if options_.hp_filter
            title = [title ' (HP filter, lambda = ' num2str(options_.hp_filter) ')'];
        end
        labels = deblank(M_.endo_names(ivar(i1),:));
        headers = char('Variables',labels);
        lh = size(labels,2)+2;
        dyntable(title,headers,labels,corr,lh,8,4);
    end
end
if options_.ar > 0 && size(stationary_vars, 1) > 0
    z=[];
    for i=1:options_.ar
        oo_.autocorr{i} = oo_.gamma_y{i+1};
        z(:,i) = diag(oo_.gamma_y{i+1}(i1,i1));
    end
    if ~options_.noprint,      
        skipline()    
        if options_.order == 2
            title='APPROXIMATED COEFFICIENTS OF AUTOCORRELATION';            
        else
            title='COEFFICIENTS OF AUTOCORRELATION';
        end
        if options_.hp_filter        
            title = [title ' (HP filter, lambda = ' num2str(options_.hp_filter) ')'];      
        end      
        labels = deblank(M_.endo_names(ivar(i1),:));      
        headers = char('Order ',int2str([1:options_.ar]'));
        lh = size(labels,2)+2;
        dyntable(title,headers,labels,z,lh,8,4);
    end  
end
