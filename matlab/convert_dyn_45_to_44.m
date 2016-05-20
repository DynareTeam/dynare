function oo_ = convert_dyn_45_to_44(M_, options_, oo_,bayestopt_)
%function oo_ = convert_dyn_45_to_44(M_, options_, oo_,bayestopt_)
% Converts oo_ from 4.5 to 4.4
%
% INPUTS
%    M_          [struct]    dynare model struct
%    options_    [struct]    dynare options struct
%    oo_         [struct]    dynare output struct
%    bayestopt_  [struct]    structure storing information about priors
% OUTPUTS
%    oo_         [struct]    dynare output struct
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2015-2016 Dynare Team
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
% along = with Dynare.  If not, see <http://www.gnu.org/licenses/>.

%% add initial conditions to Bayesian forecasts
if isfield(oo_,'PointForecast')
    var_names=fieldnames(oo_.PointForecast.HPDinf);
    moment_names=fieldnames(oo_.PointForecast);
    for moment_iter=1:length(moment_names)
        for var_iter=1:length(var_names)
            if strcmp(moment_names{moment_iter},'deciles')
                oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(:,end)*ones(M_.maximum_endo_lag,1)  oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})];
                oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(:,end)*ones(M_.maximum_endo_lag,1) oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})];                
            else
                oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(end)*ones(M_.maximum_endo_lag,1); oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})];
                oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(end)*ones(M_.maximum_endo_lag,1); oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})];
            end
        end
    end
end

%% change HPD-fields back to row vectors
if isfield(oo_,'PointForecast') && isfield(oo_.PointForecast,'HPDinf')
    names=fieldnames(oo_.PointForecast.HPDinf);
    for ii=1:length(names)
        oo_.PointForecast.HPDinf.(names{ii})=oo_.PointForecast.HPDinf.(names{ii})';
        oo_.PointForecast.HPDsup.(names{ii})=oo_.PointForecast.HPDsup.(names{ii})';
    end
end

if isfield(oo_,'MeanForecast') && isfield(oo_.MeanForecast,'HPDinf')
    names=fieldnames(oo_.MeanForecast.HPDinf);
    for ii=1:length(names)
        oo_.MeanForecast.HPDinf.(names{ii})=oo_.MeanForecast.HPDinf.(names{ii})';
        oo_.MeanForecast.HPDsup.(names{ii})=oo_.MeanForecast.HPDsup.(names{ii})';
    end
end

if isfield(oo_,'UpdatedVariables') && isfield(oo_.UpdatedVariables,'HPDinf')
    names=fieldnames(oo_.UpdatedVariables.HPDinf);
    for ii=1:length(names)
        oo_.UpdatedVariables.HPDinf.(names{ii})=oo_.UpdatedVariables.HPDinf.(names{ii})';
        oo_.UpdatedVariables.HPDsup.(names{ii})=oo_.UpdatedVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_,'SmoothedVariables') && isfield(oo_.SmoothedVariables,'HPDinf')
    names=fieldnames(oo_.SmoothedVariables.HPDinf);
    for ii=1:length(names)
        oo_.SmoothedVariables.HPDinf.(names{ii})=oo_.SmoothedVariables.HPDinf.(names{ii})';
        oo_.SmoothedVariables.HPDsup.(names{ii})=oo_.SmoothedVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_,'FilteredVariables') && isfield(oo_.FilteredVariables,'HPDinf')
    names=fieldnames(oo_.FilteredVariables.HPDinf);
    for ii=1:length(names)
        oo_.FilteredVariables.HPDinf.(names{ii})=oo_.FilteredVariables.HPDinf.(names{ii})';
        oo_.FilteredVariables.HPDsup.(names{ii})=oo_.FilteredVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_,'SmoothedShocks') && isfield(oo_.SmoothedShocks,'HPDinf')
    names=fieldnames(oo_.SmoothedShocks.HPDinf);
    for ii=1:length(names)
        oo_.SmoothedShocks.HPDinf.(names{ii})=oo_.SmoothedShocks.HPDinf.(names{ii})';
        oo_.SmoothedShocks.HPDsup.(names{ii})=oo_.SmoothedShocks.HPDsup.(names{ii})';
    end
end

%% subtract mean from classical Updated variables
if isfield(oo_,'UpdatedVariables')
    names=fieldnames(oo_.UpdatedVariables);
    for ii=1:length(names)
        %make sure Bayesian fields are not affected
        if ~strcmp(names{ii},'Mean') && ~strcmp(names{ii},'Median') && ~strcmp(names{ii},'deciles') ...
                && ~strcmp(names{ii},'Var') && ~strcmp(names{ii},'HPDinf') && ~strcmp(names{ii},'HPDsup') 
            current_var_index=find(strmatch(names{ii},deblank(M_.endo_names),'exact'));
            if  options_.loglinear == 1 %logged steady state must be used
                constant_current_variable=log(oo_.dr.ys(current_var_index));
            elseif options_.loglinear == 0 %unlogged steady state must be used
                constant_current_variable=oo_.dr.ys(current_var_index);
            end
            oo_.UpdatedVariables.(names{ii})=oo_.UpdatedVariables.(names{ii})-constant_current_variable;
            if isfield(oo_.Smoother,'Trend') && isfield(oo_.Smoother.Trend,names{ii})
                oo_.UpdatedVariables.(names{ii})=oo_.UpdatedVariables.(names{ii})-oo_.Smoother.Trend.(names{ii});
            end
        end
    end
end

%% padd classical filtered variables with redundant zeros and subtract mean
if isfield(oo_,'FilteredVariables')
    names=fieldnames(oo_.FilteredVariables);
    for ii=1:length(names)
        %make sure Bayesian fields are not affected
        if ~strcmp(names{ii},'Mean') && ~strcmp(names{ii},'Median') && ~strcmp(names{ii},'deciles') ...
                && ~strcmp(names{ii},'Var') && ~strcmp(names{ii},'HPDinf') && ~strcmp(names{ii},'HPDsup') 
            current_var_index=find(strmatch(names{ii},deblank(M_.endo_names),'exact'));
            if  options_.loglinear == 1 %logged steady state must be used
                constant_current_variable=log(oo_.dr.ys(current_var_index));
            elseif options_.loglinear == 0 %unlogged steady state must be used
                constant_current_variable=oo_.dr.ys(current_var_index);
            end
            oo_.FilteredVariables.(names{ii})=oo_.FilteredVariables.(names{ii})-constant_current_variable;
            if isfield(oo_.Smoother,'Trend') && isfield(oo_.Smoother.Trend,names{ii})
                oo_.FilteredVariables.(names{ii})=oo_.FilteredVariables.(names{ii})-oo_.Smoother.Trend.(names{ii});
            end
            oo_.FilteredVariables.(names{ii})=[0; oo_.FilteredVariables.(names{ii}); zeros(options_.nk-1,1)];
        end
    end
end

%% resort fields that are in declaration order to decision rule order (previous undocumented behavior)
if ~isempty(options_.nk) && options_.nk ~= 0 && ~isempty(bayestopt_)
    if ~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file)) %no Bayesian estimation
        positions_in_decision_order=oo_.dr.inv_order_var(bayestopt_.smoother_var_list(bayestopt_.smoother_saved_var_list));
        if ~(options_.selected_variables_only && ~(options_.forecast > 0)) %happens only when selected_variables_only is not used
            oo_.FilteredVariablesKStepAhead(:,positions_in_decision_order,:)=oo_.FilteredVariablesKStepAhead;
            if ~isempty(PK) %get K-step ahead variances
                oo_.FilteredVariablesKStepAheadVariances(:,positions_in_decision_order,positions_in_decision_order,:)=oo_.FilteredVariablesKStepAheadVariances;
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition(:,positions_in_decision_order,:,:)=oo_.FilteredVariablesShockDecomposition;
            end
        else
            fprintf('\nconvert_dyn_45_to_44:: Due to a bug in Dynare 4.4.3 with the selected_variables_only option, the previous behavior\n')
            fprintf('convert_dyn_45_to_44:: cannot be restored for FilteredVariablesKStepAhead, FilteredVariablesKStepAheadVariances, and\n')
            fprintf('convert_dyn_45_to_44:: FilteredVariablesShockDecomposition\n')
        end
    end
end

if options_.filter_covariance
    oo_.Smoother.Variance(oo_.dr.inv_order_var,oo_.dr.inv_order_var,:)=oo_.Smoother.Variance;
end


%% set old field posterior_std and remove new field posterior_std_at_mode
if isfield(oo_,'posterior_std_at_mode')
    oo_.posterior_std=oo_.posterior_std_at_mode;
    oo_=rmfield(oo_,'posterior_std_at_mode');
end


%Deal with OSR
if ~isempty(M_.osr.variable_weights)
   evalin('base','optim_weights_=M_.osr.variable_weights')
end
if ~isempty(M_.osr.variable_indices)
   evalin('base','obj_var_=M_.osr.variable_indices')
end
if ~isempty(M_.osr.param_names)
   evalin('base','osr_params_=char(M_.osr.param_names)')
end