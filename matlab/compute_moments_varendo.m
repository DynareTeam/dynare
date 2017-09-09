function oo_ = compute_moments_varendo(type,options_,M_,oo_,var_list_)
% Computes the second order moments (autocorrelation function, covariance
% matrix and variance decomposition) distributions for all the endogenous variables selected in
% var_list_. The results are saved in oo_
%
% INPUTS:
%   type            [string]       'posterior' or 'prior'
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   oo_             [structure]    Dynare structure (results).
%   var_list_       [string]       Array of string with endogenous variable names.
%
% OUTPUTS
%   oo_             [structure]    Dynare structure (results).
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2008-2017 Dynare Team
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


fprintf('Estimation::compute_moments_varendo: I''m computing endogenous moments (this may take a while)... \n');

if strcmpi(type,'posterior')
    posterior = 1;
    if nargin==4
        var_list_ = char(options_.varobs);
    end
elseif strcmpi(type,'prior')
    posterior = 0;
    if nargin==4
        var_list_ = options_.prior_analysis_endo_var_list;
        if isempty(var_list_)
            options_.prior_analysis_var_list = char(options_.varobs);
        end
    end
else
    error('compute_moments_varendo:: Unknown type!')
end

NumberOfEndogenousVariables = rows(var_list_);
NumberOfExogenousVariables = M_.exo_nbr;
NumberOfLags = options_.ar;
NoDecomposition = options_.nodecomposition;
if isfield(options_,'conditional_variance_decomposition')
    Steps = options_.conditional_variance_decomposition;
else
    Steps = 0;
end

if options_.TeX
    var_list_tex='';
    for var_iter=1:size(var_list_,1)
        var_list_tex=strvcat(var_list_tex,M_.endo_names_tex(strmatch(var_list_(var_iter,:),M_.endo_names,'exact'),:));
    end
end

% COVARIANCE MATRIX.
if posterior
    for i=1:NumberOfEndogenousVariables
        for j=i:NumberOfEndogenousVariables
            oo_ = posterior_analysis('variance',var_list_(i,:),var_list_(j,:),[],options_,M_,oo_);
        end
    end
else
    for i=1:NumberOfEndogenousVariables
        for j=i:NumberOfEndogenousVariables
            oo_ = prior_analysis('variance',var_list_(i,:),var_list_(j,:),[],options_,M_,oo_);
        end
    end
end
% CORRELATION FUNCTION.
if posterior
    for h=NumberOfLags:-1:1
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfEndogenousVariables
                oo_ = posterior_analysis('correlation',var_list_(i,:),var_list_(j,:),h,options_,M_,oo_);
            end
        end
    end
else
    for h=NumberOfLags:-1:1
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfEndogenousVariables
                oo_ = prior_analysis('correlation',var_list_(i,:),var_list_(j,:),h,options_,M_,oo_);
            end
        end
    end
end
% VARIANCE DECOMPOSITION.
if M_.exo_nbr > 1
    if ~NoDecomposition
        temp=NaN(NumberOfEndogenousVariables,NumberOfExogenousVariables);
        if posterior
            for i=1:NumberOfEndogenousVariables
                for j=1:NumberOfExogenousVariables
                    oo_ = posterior_analysis('decomposition',var_list_(i,:),M_.exo_names(j,:),[],options_,M_,oo_);
                    temp(i,j)=oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.(deblank(var_list_(i,:))).(deblank(M_.exo_names(j,:)));
                end
            end
            title='Posterior mean variance decomposition (in percent)';
            save_name_string='dsge_post_mean_var_decomp_uncond';
        else
            for i=1:NumberOfEndogenousVariables
                for j=1:NumberOfExogenousVariables
                    oo_ = prior_analysis('decomposition',var_list_(i,:),M_.exo_names(j,:),[],options_,M_,oo_);
                    temp(i,j)=oo_.PriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.(deblank(var_list_(i,:))).(deblank(M_.exo_names(j,:)));
                end
            end
            title='Prior mean variance decomposition (in percent)';
            save_name_string='dsge_prior_mean_var_decomp_uncond';
        end
        title=add_filter_subtitle(title,options_);
        headers = M_.exo_names;
        headers(M_.exo_names_orig_ord,:) = headers;
        headers = char(' ',headers);
        lh = size(deblank(var_list_),2)+2;
        dyntable(options_,title,headers,deblank(var_list_),100* ...
                 temp,lh,8,2);
        if options_.TeX
            headers=M_.exo_names_tex;
            headers = char(' ',headers);
            labels = deblank(var_list_tex);
            lh = size(labels,2)+2;
            dyn_latex_table(M_,options_,title,save_name_string,headers,labels,100*temp,lh,8,2);
        end
        skipline();
    end
    skipline();
    if ~all(M_.H==0)
        [observable_name_requested_vars,varlist_pos]=intersect(var_list_,options_.varobs,'stable');
        if ~isempty(observable_name_requested_vars)
            NumberOfObservedEndogenousVariables=length(observable_name_requested_vars);
            temp=NaN(NumberOfObservedEndogenousVariables,NumberOfExogenousVariables+1);
            if posterior
                for i=1:NumberOfObservedEndogenousVariables
                    for j=1:NumberOfExogenousVariables
                        temp(i,j,:)=oo_.PosteriorTheoreticalMoments.dsge.VarianceDecompositionME.Mean.(deblank(observable_name_requested_vars{i,1})).(deblank(M_.exo_names(j,:)));
                    end
                    endo_index_varlist=strmatch(deblank(observable_name_requested_vars{i,1}),var_list_,'exact');
                    oo_ = posterior_analysis('decomposition',var_list_(endo_index_varlist,:),'ME',[],options_,M_,oo_);
                    temp(i,j+1,:)=oo_.PosteriorTheoreticalMoments.dsge.VarianceDecompositionME.Mean.(deblank(observable_name_requested_vars{i,1})).('ME');
                end
                title='Posterior mean variance decomposition (in percent) with measurement error';
                save_name_string='dsge_post_mean_var_decomp_uncond_ME';
            else
                for i=1:NumberOfObservedEndogenousVariables
                    for j=1:NumberOfExogenousVariables
                        temp(i,j,:)=oo_.PriorTheoreticalMoments.dsge.VarianceDecompositionME.Mean.(deblank(observable_name_requested_vars(i,:))).(deblank(M_.exo_names(j,:)));
                    end
                    endo_index_varlist=strmatch(deblank(observable_name_requested_vars{i,1}),var_list_,'exact');
                    oo_ = prior_analysis('decomposition',var_list_(endo_index_varlist,:),'ME',[],options_,M_,oo_);
                    temp(i,j+1,:)=oo_.PriorTheoreticalMoments.dsge.VarianceDecompositionME.Mean.(deblank(observable_name_requested_vars{i,1})).('ME');
                end
                title='Prior mean variance decomposition (in percent) with measurement error';
                save_name_string='dsge_prior_mean_var_decomp_uncond_ME';
            end
            title=add_filter_subtitle(title,options_);
            headers = M_.exo_names;
            headers(M_.exo_names_orig_ord,:) = headers;
            headers = char(' ',headers,'ME');
            lh = size(deblank(var_list_),2)+2;
            dyntable(options_,title,headers,deblank(char(observable_name_requested_vars)),100* ...
                temp,lh,8,2);
            if options_.TeX
                headers=M_.exo_names_tex;
                headers = char(' ',headers,'ME');
                labels = deblank(var_list_tex(varlist_pos,:));
                lh = size(labels,2)+2;
                dyn_latex_table(M_,options_,title,save_name_string,headers,labels,100*temp,lh,8,2);
            end
            skipline();
        end
    end
% CONDITIONAL VARIANCE DECOMPOSITION.
    if Steps
        temp=NaN(NumberOfEndogenousVariables,NumberOfExogenousVariables,length(Steps));
        if posterior
            for i=1:NumberOfEndogenousVariables
                for j=1:NumberOfExogenousVariables
                    oo_ = posterior_analysis('conditional decomposition',i,M_.exo_names(j,:),Steps,options_,M_,oo_);
                    temp(i,j,:)=oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean.(deblank(var_list_(i,:))).(deblank(M_.exo_names(j,:)));
                end
            end
            title='Posterior mean conditional variance decomposition (in percent)';
            save_name_string='dsge_post_mean_var_decomp_cond_h';
        else
            for i=1:NumberOfEndogenousVariables
                for j=1:NumberOfExogenousVariables
                    oo_ = prior_analysis('conditional decomposition',var_list_(i,:),M_.exo_names(j,:),Steps,options_,M_,oo_);
                    temp(i,j,:)=oo_.PriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean.(deblank(var_list_(i,:))).(deblank(M_.exo_names(j,:)));
                end
            end
            title='Prior mean conditional variance decomposition (in percent)';
            save_name_string='dsge_prior_mean_var_decomp_cond_h';
        end
        for step_iter=1:length(Steps)
            title_print=[title, ' Period ' int2str(Steps(step_iter))];
            headers = M_.exo_names;
            headers(M_.exo_names_orig_ord,:) = headers;
            headers = char(' ',headers);
            lh = size(deblank(var_list_),2)+2;
            dyntable(options_,title_print,headers,deblank(var_list_),100* ...
                     temp(:,:,step_iter),lh,8,2);
            if options_.TeX
                headers=M_.exo_names_tex;
                headers = char(' ',headers);
                labels = deblank(var_list_tex);
                lh = size(labels,2)+2;
                dyn_latex_table(M_,options_,title_print,[save_name_string,int2str(Steps(step_iter))],headers,labels,100*temp(:,:,step_iter),lh,8,2);
            end
        end
        skipline();
        if ~all(M_.H==0)
            if ~isempty(observable_name_requested_vars)
                NumberOfObservedEndogenousVariables=length(observable_name_requested_vars);
                temp=NaN(NumberOfObservedEndogenousVariables,NumberOfExogenousVariables+1,length(Steps));
                if posterior
                    for i=1:NumberOfObservedEndogenousVariables
                        for j=1:NumberOfExogenousVariables
                            temp(i,j,:)=oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecompositionME.Mean.(deblank(observable_name_requested_vars{i,1})).(deblank(M_.exo_names(j,:)));
                        end
                        endo_index_varlist=strmatch(deblank(observable_name_requested_vars{i,1}),var_list_,'exact');
                        oo_ = posterior_analysis('conditional decomposition',endo_index_varlist,'ME',Steps,options_,M_,oo_);
                        temp(i,j+1,:)=oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecompositionME.Mean.(deblank(observable_name_requested_vars{i,1})).('ME');
                    end                    
                    title='Posterior mean conditional variance decomposition (in percent) with measurement error';
                    save_name_string='dsge_post_mean_var_decomp_ME_cond_h';
                else
                   for i=1:NumberOfObservedEndogenousVariables
                        for j=1:NumberOfExogenousVariables
                            temp(i,j,:)=oo_.PriorTheoreticalMoments.dsge.ConditionalVarianceDecompositionME.Mean.(deblank(observable_name_requested_vars(i,:))).(deblank(M_.exo_names(j,:)));
                        end
                        endo_index_varlist=strmatch(deblank(observable_name_requested_vars{i,1}),var_list_,'exact');
                        oo_ = prior_analysis('conditional decomposition',endo_index_varlist,'ME',Steps,options_,M_,oo_);
                        temp(i,j+1,:)=oo_.PriorTheoreticalMoments.dsge.ConditionalVarianceDecompositionME.Mean.(deblank(observable_name_requested_vars{i,1})).('ME');
                    end                    
                    title='Prior mean conditional variance decomposition (in percent) with measurement error';
                    save_name_string='dsge_prior_mean_var_decomp_ME_cond_h';
                end
                for step_iter=1:length(Steps)
                    title_print=[title, ' Period ' int2str(Steps(step_iter))];
                    headers = M_.exo_names;
                    headers(M_.exo_names_orig_ord,:) = headers;
                    headers = char(' ',headers,'ME');
                    lh = size(deblank(var_list_),2)+2;
                    dyntable(options_,title_print,headers,deblank(char(observable_name_requested_vars)),100* ...
                        temp(:,:,step_iter),lh,8,2);
                    if options_.TeX
                        headers=M_.exo_names_tex;
                        headers = char(' ',headers,'ME');
                        labels = deblank(var_list_tex(varlist_pos,:));
                        lh = size(labels,2)+2;
                        dyn_latex_table(M_,options_,title_print,[save_name_string,int2str(Steps(step_iter))],headers,labels,100*temp(:,:,step_iter),lh,8,2);
                    end
                end
                skipline();
            end
        end        
    end
end

fprintf(' Done!\n');
