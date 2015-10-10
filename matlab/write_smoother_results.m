function [oo_, yf]=write_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend)
% oo_=write_smoother_results(M_,oo_,options_,bayestopt_,dataset_,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend)
% Writes the smoother results into respective fields in oo_
% 
% Inputs:
%   M_              [structure]     storing the model information
%   oo_             [structure]     storing the results
%   options_        [structure]     storing the options
%   bayestopt_      [structure]     storing information about priors
%   dataset_        [structure]     storing the dataset
%   atT             [double]    (m*T) matrix, smoothed endogenous variables (a_{t|T})
%   innov           [double]    (r*T) matrix, smoothed structural shocks (r>n is the umber of shocks).
%   measurement_error [double]  (n*T) matrix, smoothed measurement errors.
%   updated_variables [double]  (m*T) matrix, updated (endogenous) variables (a_{t|T})
%   ys              [double]    (m*1) vector specifying the steady state level of each endogenous variable.
%   trend_coeff     [double]    (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
%   aK              [double]    (K,n,T+K) array, k (k=1,...,K) steps ahead filtered (endogenous) variables.
%   P               [3D array]  of one-step ahead forecast error variance
%                   matrices
%   PK              [4D array]  of k-step ahead forecast error variance
%                               matrices (meaningless for periods 1:d)
%   decomp
%   Trend           [double]    [nvarobs*T] matrix of trends in observables
%
% Outputs:
%   oo_             [structure] storing the results
%   yf              [double]    (nvarobs*T) matrix storing the smoothed observed variables   
% 
% Notes: first all smoothed variables are saved without trend and constant.
%       Then trend and constant are added for the observed variables.
%
% Copyright (C) 2014 Dynare Team
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
oo_.Smoother.TrendCoeffs = trend_coeff; %are in order of options_.varobs

if ~isempty(Trend)
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

%get indicees of smoothed variables
i_endo = bayestopt_.smoother_saved_var_list;

if ~isempty(options_.nk) && options_.nk ~= 0 && (~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file)))
    %write deviations from steady state, add constant for observables later
    oo_.FilteredVariablesKStepAhead = aK(options_.filter_step_ahead,i_endo,:);    
    if ~isempty(PK) %get K-step ahead variances
        oo_.FilteredVariablesKStepAheadVariances = ...
            PK(options_.filter_step_ahead,i_endo,i_endo,:);
    end
    if ~isempty(decomp) %get decomposition
        oo_.FilteredVariablesShockDecomposition = ...
            decomp(options_.filter_step_ahead,i_endo,:,:);
    end
end

for i=bayestopt_.smoother_saved_var_list'
    i1 = oo_.dr.order_var(bayestopt_.smoother_var_list(i)); %get indices of smoothed variables in name vector
    %% Compute constant
    if  options_.loglinear == 1 %logged steady state must be used
        constant_current_variable=repmat(log(ys(i1)),gend,1);
    elseif options_.loglinear == 0 %unlogged steady state must be used
        constant_current_variable=repmat((ys(i1)),gend,1);
    end
    oo_.SmoothedVariables.(deblank(M_.endo_names(i1,:)))=atT(i,:)'+constant_current_variable;
    if ~isempty(options_.nk) && options_.nk > 0 && ~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file))
        oo_.FilteredVariables.(deblank(M_.endo_names(i1,:)))=squeeze(aK(1,i,2:end-(options_.nk-1)));
    end
    oo_.UpdatedVariables.(deblank(M_.endo_names(i1,:)))=updated_variables(i,:)'+constant_current_variable;
end
    
%% Add trend and constant for observed variables
for pos_iter=1:length(bayestopt_.mf)
    oo_.Smoother.Constant.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=constant_part(pos_iter,:);
    oo_.SmoothedVariables.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=yf(pos_iter,:)';   
    if ~isempty(options_.nk) && options_.nk > 0 && ~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file))
        %filtered variable E_t(y_t+1) requires to shift trend by 1 period
        oo_.FilteredVariables.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=...
            squeeze(aK(1,bayestopt_.mf(pos_iter),2:end-(options_.nk-1)))...
            +trend_constant_observables_filtered.filter_ahead_1(pos_iter,:)';
        for filter_iter=1:length(options_.filter_step_ahead)
            filter_step=options_.filter_step_ahead(filter_iter);
            oo_.FilteredVariablesKStepAhead(filter_iter,bayestopt_.mf(pos_iter),1+filter_step:end-(max(options_.filter_step_ahead)-filter_step)) = ...
                squeeze(aK(filter_step,bayestopt_.mf(pos_iter),1+filter_step:end-(max(options_.filter_step_ahead)-filter_step)))...
                +trend_constant_observables_filtered.(['filter_ahead_' num2str(filter_step)])(pos_iter,:)';    
        end
    end
    %updated variables are E_t(y_t) so no trend shift is required
    oo_.UpdatedVariables.(deblank(M_.endo_names(bayestopt_.mfys(pos_iter),:)))=...
        updated_variables(bayestopt_.mf(pos_iter),:)'+trend_constant_observables(pos_iter,:)';
end

%% get smoothed shocks
for i=1:M_.exo_nbr
    oo_.SmoothedShocks.(deblank(M_.exo_names(i,:)))=innov(i,:)';
end

%%  Smoothed measurement errors
if ~isequal(M_.H,0)
%     measurement_error_indices=find(diag(M_.H)~=0);
    for meas_error_iter=1:length(options_.varobs)
       oo_.SmoothedMeasurementErrors.(options_.varobs{meas_error_iter})= measurement_error(meas_error_iter,:)';
    end
end