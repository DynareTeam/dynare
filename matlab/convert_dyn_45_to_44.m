function []=convert_dyn_45_to_44
% function []=convert_dyn_45_to_44
% This function converts the oo_-structure fields that have been changed in Dynare 4.5. 
% following https://github.com/DynareTeam/dynare/pull/771 to the old format
% of Dynare 4.4

global M_ oo_ options_

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
if isfield(oo_.PointForecast,'HPDinf')
    names=fieldnames(oo_.PointForecast.HPDinf);
    for ii=1:length(names)
        oo_.PointForecast.HPDinf.(names{ii})=oo_.PointForecast.HPDinf.(names{ii})';
        oo_.PointForecast.HPDsup.(names{ii})=oo_.PointForecast.HPDsup.(names{ii})';
    end
end

if isfield(oo_.MeanForecast,'HPDinf')
    names=fieldnames(oo_.MeanForecast.HPDinf);
    for ii=1:length(names)
        oo_.MeanForecast.HPDinf.(names{ii})=oo_.MeanForecast.HPDinf.(names{ii})';
        oo_.MeanForecast.HPDsup.(names{ii})=oo_.MeanForecast.HPDsup.(names{ii})';
    end
end

if isfield(oo_.UpdatedVariables,'HPDinf')
    names=fieldnames(oo_.UpdatedVariables.HPDinf);
    for ii=1:length(names)
        oo_.UpdatedVariables.HPDinf.(names{ii})=oo_.UpdatedVariables.HPDinf.(names{ii})';
        oo_.UpdatedVariables.HPDsup.(names{ii})=oo_.UpdatedVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_.SmoothedVariables,'HPDinf')
    names=fieldnames(oo_.SmoothedVariables.HPDinf);
    for ii=1:length(names)
        oo_.SmoothedVariables.HPDinf.(names{ii})=oo_.SmoothedVariables.HPDinf.(names{ii})';
        oo_.SmoothedVariables.HPDsup.(names{ii})=oo_.SmoothedVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_.FilteredVariables,'HPDinf')
    names=fieldnames(oo_.FilteredVariables.HPDinf);
    for ii=1:length(names)
        oo_.FilteredVariables.HPDinf.(names{ii})=oo_.FilteredVariables.HPDinf.(names{ii})';
        oo_.FilteredVariables.HPDsup.(names{ii})=oo_.FilteredVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_.SmoothedShocks,'HPDinf')
    names=fieldnames(oo_.SmoothedShocks.HPDinf);
    for ii=1:length(names)
        oo_.SmoothedShocks.HPDinf.(names{ii})=oo_.SmoothedShocks.HPDinf.(names{ii})';
        oo_.SmoothedShocks.HPDsup.(names{ii})=oo_.SmoothedShocks.HPDsup.(names{ii})';
    end
end

%% padd classical filtered variables with redundant zeros
if isfield(oo_,'FilteredVariables')
    names=fieldnames(oo_.FilteredVariables);
    for ii=1:length(names)
        %make sure Bayesian fields are not affect
        if ~strcmp(names{ii},'Mean') && ~strcmp(names{ii},'Median') && ~strcmp(names{ii},'deciles') ...
                && ~strcmp(names{ii},'Var') && ~strcmp(names{ii},'HPDinf') && ~strcmp(names{ii},'HPDsup') 
            oo_.FilteredVariables.(names{ii})=[0; oo_.FilteredVariables.(names{ii}); zeros(options_.nk-1,1)];
        end
    end
end
