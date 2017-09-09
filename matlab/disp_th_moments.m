function oo_=disp_th_moments(dr,var_list,M_,options_,oo_)
% Display theoretical moments of variables

% Copyright (C) 2001-2017 Dynare Team
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

nodecomposition = options_.nodecomposition;
if options_.one_sided_hp_filter
    error(['disp_th_moments:: theoretical moments incompatible with one-sided HP filter. Use simulated moments instead'])
end
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

[oo_.gamma_y,stationary_vars] = th_autocovariances(dr,ivar,M_,options_, nodecomposition);
m = dr.ys(ivar);
non_stationary_vars = setdiff(1:length(ivar),stationary_vars);
m(non_stationary_vars) = NaN;

i1 = find(abs(diag(oo_.gamma_y{1})) > 1e-12);
s2 = diag(oo_.gamma_y{1});
sd = sqrt(s2);
if options_.order == 2 && ~M_.hessian_eq_zero
    m = m+oo_.gamma_y{options_.ar+3};
end

z = [ m sd s2 ];
oo_.mean = m;
oo_.var = oo_.gamma_y{1};

ME_present=0;
if ~all(M_.H==0)
    [observable_pos_requested_vars,index_subset,index_observables]=intersect(ivar,options_.varobs_id,'stable');
    if ~isempty(observable_pos_requested_vars)
        ME_present=1;
    end
end
                
if size(stationary_vars, 1) > 0
    if ~nodecomposition
        oo_.variance_decomposition=100*oo_.gamma_y{options_.ar+2};
        if ME_present
            ME_Variance=diag(M_.H);
            oo_.variance_decomposition_ME=oo_.variance_decomposition(index_subset,:).*repmat(diag(oo_.var(index_subset,index_subset))./(diag(oo_.var(index_subset,index_subset))+ME_Variance(index_observables)),1,M_.exo_nbr);
            oo_.variance_decomposition_ME(:,end+1)=100-sum(oo_.variance_decomposition_ME,2);
        end
    end
    if ~options_.noprint %options_.nomoments == 0
        if options_.order == 2
            title='APPROXIMATED THEORETICAL MOMENTS';
        else
            title='THEORETICAL MOMENTS';
        end
        title=add_filter_subtitle(title,options_);
        headers=char('VARIABLE','MEAN','STD. DEV.','VARIANCE');
        labels = deblank(M_.endo_names(ivar,:));
        lh = size(labels,2)+2;
        dyntable(options_,title,headers,labels,z,lh,11,4);
        if options_.TeX
            labels = deblank(M_.endo_names_tex(ivar,:));
            lh = size(labels,2)+2;
            dyn_latex_table(M_,options_,title,'th_moments',headers,labels,z,lh,11,4);
        end

        if M_.exo_nbr > 1 && ~nodecomposition
            skipline()
            if options_.order == 2
                title='APPROXIMATED VARIANCE DECOMPOSITION (in percent)';
            else
                title='VARIANCE DECOMPOSITION (in percent)';
            end
            title=add_filter_subtitle(title,options_);
            headers = M_.exo_names;
            headers(M_.exo_names_orig_ord,:) = headers;
            headers = char(' ',headers);
            lh = size(deblank(M_.endo_names(ivar(stationary_vars),:)),2)+2;
            dyntable(options_,title,headers,deblank(M_.endo_names(ivar(stationary_vars), ...
                                                              :)),100* ...
                     oo_.gamma_y{options_.ar+2}(stationary_vars,:),lh,8,2);
            if ME_present
                [stationary_observables,pos_index_subset]=intersect(index_subset,stationary_vars,'stable');
                headers_ME=char(headers,'ME');
                dyntable(options_,[title,' WITH MEASUREMENT ERROR'],headers_ME,deblank(M_.endo_names(ivar(stationary_observables), ...
                                                  :)),oo_.variance_decomposition_ME(pos_index_subset,:),lh,8,2);
            end
            if options_.TeX
                headers=M_.exo_names_tex;
                headers = char(' ',headers);
                labels = deblank(M_.endo_names_tex(ivar(stationary_vars),:));
                lh = size(labels,2)+2;
                dyn_latex_table(M_,options_,title,'th_var_decomp_uncond',headers,labels,100*oo_.gamma_y{options_.ar+2}(stationary_vars,:),lh,8,2);
                if ME_present
                    headers_ME=char(headers,'ME');
                    dyn_latex_table(M_,options_,[title,' WITH MEASUREMENT ERROR'],'th_var_decomp_uncond_ME',headers_ME,labels,oo_.variance_decomposition_ME(pos_index_subset,:),lh,8,2);
                end
            end
        end
    end

    conditional_variance_steps = options_.conditional_variance_decomposition;
    if length(conditional_variance_steps)
        StateSpaceModel.number_of_state_equations = M_.endo_nbr;
        StateSpaceModel.number_of_state_innovations = M_.exo_nbr;
        StateSpaceModel.sigma_e_is_diagonal = M_.sigma_e_is_diagonal;
        [StateSpaceModel.transition_matrix,StateSpaceModel.impulse_matrix] = kalman_transition_matrix(dr,(1:M_.endo_nbr)',M_.nstatic+(1:M_.nspred)',M_.exo_nbr);
        StateSpaceModel.state_innovations_covariance_matrix = M_.Sigma_e;
        StateSpaceModel.order_var = dr.order_var;
        StateSpaceModel.measurement_error=M_.H;
        StateSpaceModel.observable_pos=options_.varobs_id;
        [oo_.conditional_variance_decomposition, oo_.conditional_variance_decomposition_ME]= conditional_variance_decomposition(StateSpaceModel,conditional_variance_steps,ivar);

        if options_.noprint == 0
            display_conditional_variance_decomposition(oo_.conditional_variance_decomposition,conditional_variance_steps,...
                                                       ivar,M_,options_);
            if ME_present
                display_conditional_variance_decomposition(oo_.conditional_variance_decomposition_ME,conditional_variance_steps,...
                                                       observable_pos_requested_vars,M_,options_);               
            end                                                       
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
    corr=NaN(size(oo_.gamma_y{1}));
    corr(i1,i1) = oo_.gamma_y{1}(i1,i1)./(sd(i1)*sd(i1)');
    if options_.contemporaneous_correlation 
        oo_.contemporaneous_correlation = corr;
    end
    if ~options_.noprint
        skipline()
        if options_.order == 2
            title='APPROXIMATED MATRIX OF CORRELATIONS';
        else
            title='MATRIX OF CORRELATIONS';
        end
        title=add_filter_subtitle(title,options_);
        labels = deblank(M_.endo_names(ivar(i1),:));
        headers = char('Variables',labels);
        lh = size(labels,2)+2;
        dyntable(options_,title,headers,labels,corr(i1,i1),lh,8,4);
        if options_.TeX
            labels = deblank(M_.endo_names_tex(ivar(i1),:));
            headers=char('Variables',labels);
            lh = size(labels,2)+2;
            dyn_latex_table(M_,options_,title,'th_corr_matrix',headers,labels,corr(i1,i1),lh,8,4);
        end
    end
end
if options_.ar > 0 && size(stationary_vars, 1) > 0
    z=[];
    for i=1:options_.ar
        oo_.autocorr{i} = oo_.gamma_y{i+1};
        z(:,i) = diag(oo_.gamma_y{i+1}(i1,i1));
    end
    if ~options_.noprint
        skipline()
        if options_.order == 2
            title='APPROXIMATED COEFFICIENTS OF AUTOCORRELATION';
        else
            title='COEFFICIENTS OF AUTOCORRELATION';
        end
        title=add_filter_subtitle(title,options_);
        labels = deblank(M_.endo_names(ivar(i1),:));
        headers = char('Order ',int2str([1:options_.ar]'));
        lh = size(labels,2)+2;
        dyntable(options_,title,headers,labels,z,lh,8,4);
        if options_.TeX
            labels = deblank(M_.endo_names_tex(ivar(i1),:));
            headers=char('Order ',int2str([1:options_.ar]'));
            lh = size(labels,2)+2;
            dyn_latex_table(M_,options_,title,'th_autocorr_matrix',headers,labels,z,lh,8,4);
        end
    end
end
