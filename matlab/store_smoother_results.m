function [oo_, yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend)
% oo_=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend)
% Writes the smoother results into respective fields in oo_
% 
% Inputs:
%   M_              [structure]     storing the model information
%   oo_             [structure]     storing the results
%   options_        [structure]     storing the options
%   bayestopt_      [structure]     storing information about priors
%   dataset_        [structure]     storing the dataset
%   atT             [double]    (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
%   innov           [double]    (r*T) matrix, smoothed structural shocks (r>n is the umber of shocks).
%   measurement_error [double]  (n*T) matrix, smoothed measurement errors.
%   updated_variables [double]  (m*T) matrix, updated (endogenous) variables (a_{t|t}) (decision-rule order)
%   ys              [double]    (m*1) vector specifying the steady state
%                                   level of each endogenous variable (declaration order)
%   trend_coeff     [double]    (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
%   aK              [double]    (K,n,T+K) array, k (k=1,...,K) steps ahead
%                                   filtered (endogenous) variables  (decision-rule order)
%   P               [3D array]  (m*m*(T+1)) array of one-step ahead forecast error variance
%                                   matrices (decision-rule order)
%   PK              [4D array]  (K*m*m*(T+K)) 4D array of k-step ahead forecast error variance
%                                   matrices (meaningless for periods 1:d) (decision-rule order)
%   decomp          [4D array]  (K*m*r*(T+K)) 4D array of shock decomposition of k-step ahead
%                                   filtered variables (decision-rule order)
%   Trend           [double]    [nvarobs*T] matrix of trends in observables
%
% Outputs:
%   oo_             [structure] storing the results:
%                   oo_.Smoother.SteadyState: Steady states (declaration order)
%                   oo_.Smoother.TrendCoeffs: trend coefficients, with zeros where no trend applies (declaration order)
%                   oo_.Smoother.Variance: one-step ahead forecast error variance (declaration order)
%                   oo_.Smoother.Constant: structure storing the constant term of the smoother
%                   oo_.Smoother.Trend: structure storing the trend term of the smoother
%                   oo_.FilteredVariablesKStepAhead: k-step ahead forecast error variance matrices (decision-rule order)
%                   oo_.FilteredVariablesShockDecomposition: shock decomposition of k-step ahead filtered variables (decision-rule order)
%                   oo_.FilteredVariables: structure storing the filtered variables
%                   oo_.UpdatedVariables: structure storing the updated variables
%                   oo_.SmoothedShocks: structure storing the smoothed shocks
%                   oo_.SmoothedMeasurementErrors: structure storing the smoothed measurement errors
%
%   yf              [double]    (nvarobs*T) matrix storing the smoothed observed variables (order of options_.varobs)  
% 
% Notes: 
%   m:  number of endogenous variables (M_.endo_nbr)
%   T:  number of Time periods (options_.nobs)
%   r:  number of strucural shocks (M_.exo_nbr)
%   n:  number of observables (length(options_.varobs))
%   K:  maximum forecast horizon (max(options_.nk))
% 
%   First all smoothed variables are saved without trend and constant.
%       Then trend and constant are added for the observed variables.
%
% Copyright (C) 2014-2016 Dynare Team
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

gend=dataset_.nobs;
if nargin<16
    Trend=zeros(options_.number_of_observed_variables,gend);
end

%% write variable steady state
oo_.Smoother.SteadyState = ys;

%% write trend coefficients and trend
oo_.Smoother.TrendCoeffs = zeros(size(ys));

oo_.Smoother.TrendCoeffs(options_.varobs_id)=trend_coeff; %are in order of options_.varobs

if ~isempty(Trend)
    for var_iter=1:M_.endo_nbr
        if isempty(strmatch(deblank(M_.endo_names(var_iter,:)),options_.varobs,'exact'))
            oo_.Smoother.Trend.(deblank(M_.endo_names(var_iter,:))) = zeros(gend,1);
        end
    end
    for var_iter=1:size(options_.varobs,2)
        oo_.Smoother.Trend.(deblank(options_.varobs{1,var_iter})) = Trend(var_iter,:)';
    end
end
%% Compute constant for observables
if options_.prefilter == 1 %as mean is taken after log transformation, no distinction is needed here
    constant_part=repmat(dataset_info.descriptive.mean',1,gend);
elseif options_.prefilter == 0 && options_.loglinear == 1 %logged steady state must be used
    constant_part=repmat(log(ys(bayestopt_.mfys)),1,gend);
elseif options_.prefilter == 0 && options_.loglinear == 0 %unlogged steady state must be used
    constant_part=repmat(ys(bayestopt_.mfys),1,gend);
end

%% get observed variables including trend and constant
trend_constant_observables=constant_part+Trend;
yf = atT(bayestopt_.mf,:)+trend_constant_observables;

if options_.nk > 0
    %filtered variable E_t(y_t+k) requires to shift trend by k periods    
    filter_steps_required=union(1,options_.filter_step_ahead); % 1 is required for standard filtered variables
    for filter_iter=1:length(filter_steps_required)
        filter_step=filter_steps_required(filter_iter);
        trend_constant_observables_filtered.(['filter_ahead_' num2str(filter_step)])=constant_part+[Trend+repmat(filter_step*trend_coeff,1,gend)];
    end
end
%% write smoother variance
if options_.filter_covariance
    oo_.Smoother.Variance = P;
end

%get indices of smoothed variables
i_endo_in_bayestopt_smoother_varlist = bayestopt_.smoother_saved_var_list;
i_endo_in_dr_matrices=bayestopt_.smoother_var_list(i_endo_in_bayestopt_smoother_varlist);
if ~isempty(options_.nk) && options_.nk ~= 0
    %% Compute constant
    i_endo_declaration_order = oo_.dr.order_var(i_endo_in_dr_matrices); %get indices of smoothed variables in name vector
    if  options_.loglinear == 1 %logged steady state must be used
        constant_all_variables=repmat(log(ys(i_endo_declaration_order))',[length(options_.filter_step_ahead),1,gend+max(options_.filter_step_ahead)]);
    elseif options_.loglinear == 0 %unlogged steady state must be used
        constant_all_variables=repmat((ys(i_endo_declaration_order))',[length(options_.filter_step_ahead),1,gend+max(options_.filter_step_ahead)]);
    end
    % add constant
    oo_.FilteredVariablesKStepAhead = aK(options_.filter_step_ahead,i_endo_in_dr_matrices,:)+constant_all_variables;
    if ~isempty(PK) %get K-step ahead variances
        oo_.FilteredVariablesKStepAheadVariances = ...
            PK(options_.filter_step_ahead,i_endo_in_dr_matrices,i_endo_in_dr_matrices,:);
    end
    if ~isempty(decomp) %get decomposition
        oo_.FilteredVariablesShockDecomposition = ...
            decomp(options_.filter_step_ahead,i_endo_in_dr_matrices,:,:);
    end
end

for i_endo_in_bayestopt_smoother_varlist=bayestopt_.smoother_saved_var_list'
    i_endo_in_dr=bayestopt_.smoother_var_list(i_endo_in_bayestopt_smoother_varlist);
    i_endo_declaration_order = oo_.dr.order_var(i_endo_in_dr); %get indices of smoothed variables in name vector
    %% Compute constant
    if  options_.loglinear == 1 %logged steady state must be used
        constant_current_variable=repmat(log(ys(i_endo_declaration_order)),gend,1);
    elseif options_.loglinear == 0 %unlogged steady state must be used
        constant_current_variable=repmat((ys(i_endo_declaration_order)),gend,1);
    end
    oo_.Smoother.Constant.(deblank(M_.endo_names(i_endo_declaration_order,:)))=constant_current_variable;
    oo_.SmoothedVariables.(deblank(M_.endo_names(i_endo_declaration_order,:)))=atT(i_endo_in_dr,:)'+constant_current_variable;
    if ~isempty(options_.nk) && options_.nk > 0 && ~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file))
        oo_.FilteredVariables.(deblank(M_.endo_names(i_endo_declaration_order,:)))=squeeze(aK(1,i_endo_in_dr,2:end-(options_.nk-1)))+constant_current_variable;
    end
    oo_.UpdatedVariables.(deblank(M_.endo_names(i_endo_declaration_order,:)))=updated_variables(i_endo_in_dr,:)'+constant_current_variable;
end
    
%% Add trend and constant for observed variables
for pos_iter=1:length(bayestopt_.mf)
    oo_.Smoother.Constant.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=constant_part(pos_iter,:)';
    if ismember(bayestopt_.mf(pos_iter),bayestopt_.smoother_var_list(bayestopt_.smoother_saved_var_list))
        oo_.SmoothedVariables.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=yf(pos_iter,:)';   
        if ~isempty(options_.nk) && options_.nk > 0 
            %filtered variable E_t(y_t+1) requires to shift trend by 1 period
            oo_.FilteredVariables.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=...
                squeeze(aK(1,bayestopt_.mf(pos_iter),2:end-(options_.nk-1)))...
                +trend_constant_observables_filtered.filter_ahead_1(pos_iter,:)';
            for filter_iter=1:length(options_.filter_step_ahead)
                filter_step=options_.filter_step_ahead(filter_iter);
                oo_.FilteredVariablesKStepAhead(filter_iter,find(i_endo_in_dr_matrices==bayestopt_.mf(pos_iter)),1+filter_step:end-(max(options_.filter_step_ahead)-filter_step)) = ...
                    squeeze(aK(filter_step,bayestopt_.mf(pos_iter),1+filter_step:end-(max(options_.filter_step_ahead)-filter_step)))...
                    +trend_constant_observables_filtered.(['filter_ahead_' num2str(filter_step)])(pos_iter,:)';    
            end
        end
        %updated variables are E_t(y_t) so no trend shift is required
        oo_.UpdatedVariables.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=...
            updated_variables(bayestopt_.mf(pos_iter),:)'+trend_constant_observables(pos_iter,:)';
    end
end

%% resort fields that are in decision rule order to declaration order
if ~isempty(options_.nk) && options_.nk ~= 0
    positions_in_declaration_order=oo_.dr.order_var(bayestopt_.smoother_var_list(bayestopt_.smoother_saved_var_list));
    if ~(options_.selected_variables_only && ~(options_.forecast > 0)) %happens only when selected_variables_only is not used
        oo_.FilteredVariablesKStepAhead(:,positions_in_declaration_order,:)=oo_.FilteredVariablesKStepAhead;
        if ~isempty(PK) %get K-step ahead variances
            oo_.FilteredVariablesKStepAheadVariances(:,positions_in_declaration_order,positions_in_declaration_order,:)=oo_.FilteredVariablesKStepAheadVariances;
        end
        if ~isempty(decomp)
            oo_.FilteredVariablesShockDecomposition(:,positions_in_declaration_order,:,:)=oo_.FilteredVariablesShockDecomposition;
        end
    else
        positions_in_declaration_order=oo_.dr.order_var(bayestopt_.smoother_var_list(bayestopt_.smoother_saved_var_list));
        [junk,sorted_index_declaration_order]=sort(positions_in_declaration_order);
        oo_.FilteredVariablesKStepAhead(:,sorted_index_declaration_order,:)=oo_.FilteredVariablesKStepAhead;
        if ~isempty(PK) %get K-step ahead variances
            oo_.FilteredVariablesKStepAheadVariances(:,sorted_index_declaration_order,sorted_index_declaration_order,:)=oo_.FilteredVariablesKStepAheadVariances;
        end
        if ~isempty(decomp)
            oo_.FilteredVariablesShockDecomposition(:,sorted_index_declaration_order,:,:)=oo_.FilteredVariablesShockDecomposition;
        end
    end
end

if options_.filter_covariance
    oo_.Smoother.Variance(oo_.dr.order_var,oo_.dr.order_var,:)=oo_.Smoother.Variance;
end

%% get smoothed shocks
for exo_iter=1:M_.exo_nbr
    oo_.SmoothedShocks.(deblank(M_.exo_names(exo_iter,:)))=innov(exo_iter,:)';
end

%%  Smoothed measurement errors
if ~isequal(M_.H,0)
%     measurement_error_indices=find(diag(M_.H)~=0);
    for meas_error_iter=1:length(options_.varobs)
       oo_.SmoothedMeasurementErrors.(options_.varobs{meas_error_iter})= measurement_error(meas_error_iter,:)';
    end
end
