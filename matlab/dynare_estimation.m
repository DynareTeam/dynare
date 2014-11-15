function dynare_estimation(var_list,dname)
% function dynare_estimation(var_list)
% runs the estimation of the model
%
% INPUTS
%   var_list:  selected endogenous variables vector
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2013 Dynare Team
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

global options_ oo_ M_ oo_recursive_

% Test if the order of approximation is nonzero (the preprocessor tests if order is non negative).
if isequal(options_.order,0)
    error('Estimation:: The order of the Taylor approximation cannot be 0!')
end

% Decide if a DSGE or DSGE-VAR has to be estimated.
if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
    options_.dsge_var = 1;
end

var_list = check_list_of_variables(options_, M_, var_list);
options_.varlist = var_list;

if isfield(options_,'nobs')
    nobs = sort(options_.nobs); 
else
    nobs = [];
end

nnobs = length(nobs);
horizon = options_.forecast;

if nargin<2 || ~exist(dname) || isempty(dname)
    dname = M_.fname;
end

M_.dname = dname;

if (isnumeric(options_.mode_compute) && options_.mode_compute && options_.analytic_derivation) ... %no user supplied function
        || (~isnumeric(options_.mode_compute) && options_.analytic_derivation) % user supplied function
    analytic_derivation0=options_.analytic_derivation;
    options_.analytic_derivation=1;
end

if options_.logged_steady_state
    oo_.dr.ys=exp(oo_.dr.ys);
    oo_.steady_state=exp(oo_.steady_state);
end


if nnobs > 1
    for i=1:nnobs
        options_.nobs = nobs(i);
        M_.dname = [dname '_' int2str(nobs(i))];
        dynare_estimation_1(var_list,M_.dname);
        if isequal(i,1)
            options_.mode_file = [M_.fname '_mode'];
        end
        if options_.recursive_estimation_restart
            for j=1:options_.recursive_estimation_restart
                dynare_estimation_1(var_list,M_.dname);
            end
        end
        oo_recursive_{nobs(i)} = oo_;
    end
else
    dynare_estimation_1(var_list,dname);
end

if options_.mode_compute && options_.analytic_derivation,
    options_.analytic_derivation=analytic_derivation0;
end

if nnobs > 1 && horizon > 0
    mh_replic = options_.mh_replic;
    rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
    gend = options_.nobs;
    data_plot_end_point=min(options_.first_obs+gend-1+horizon,size(rawdata,1)); %compute last observation that can be plotted
    rawdata = rawdata(options_.first_obs:data_plot_end_point,:);
    % Take the log of the variables if needed
    if options_.loglinear && ~options_.logdata   % and if the data are not in logs, then...
        rawdata = log(rawdata);
    end

    endo_names = M_.endo_names;
    n_varobs = length(options_.varobs);

    if isempty(var_list)
        var_list = endo_names;
        nvar    = size(endo_names,1);
        SelecVariables = transpose(1:nvar);
    else
        nvar = size(var_list,1);
        SelecVariables = [];
        for i=1:nvar
            if ~isempty(strmatch(var_list(i,:),endo_names,'exact'))
                SelecVariables = [SelecVariables;strmatch(var_list(i,:),endo_names, 'exact')];
            else
                error(['Estimation:: ' var_list(i,:) ' isn''t an endogenous variable'])
            end
        end
    end

    IdObs    = zeros(n_varobs,1);
    for j=1:n_varobs
        iobs = strmatch(options_.varobs{j},var_list,'exact');
        if ~isempty(iobs)
            IdObs(j,1) = iobs;
        end
    end

    time_offset=min(3,gend-1); %for observables, plot 3 previous periods unless data is shorter
    k = time_offset+min(nobs(end)-nobs(1)+horizon, ...
              size(rawdata,1)-nobs(1));
    data2 = rawdata(end-k+1:end,:);
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(nvar);
    m = 1;
    plot_index=0;
    OutputDirectoryName = CheckPath('graphs',M_.fname);
    for i = 1:size(var_list,1)
        if mod(i,nstar) == 1
            plot_index=plot_index+1;
            hfig = dyn_figure(options_,'Name',['Out of sample forecasts (',num2str(plot_index),')']);
            m = 1;
        end
        subplot(nr,nc,m)
        hold on
        if any(i==IdObs)
            k2 = find(i==IdObs);
            offsetx = 3;
            plot(nobs(1)-offsetx+1:nobs(1)-offsetx+k,data2(end-k+1:end,k2)','-k','linewidth',2);
        else
            offsetx = 0;
        end
        vname = deblank(var_list(i,:));
        for j=1:nnobs
            if mh_replic > 0
                oo_.RecursiveForecast.Mean.(vname)(j,:) = ...
                      oo_recursive_{nobs(j)}.MeanForecast.Mean.(vname);
                oo_.RecursiveForecast.HPDinf.(vname)(j,:) = ...
                      oo_recursive_{nobs(j)}.MeanForecast.HPDinf.(vname);
                oo_.RecursiveForecast.HPDsup.(vname)(j,:) = ...
                      oo_recursive_{nobs(j)}.MeanForecast.HPDsup.(vname);
                oo_.RecursiveForecast.HPDTotalinf.(vname)(j,:) = ...
                      oo_recursive_{nobs(j)}.PointForecast.HPDinf.(vname);
                oo_.RecursiveForecast.HPDTotalsup.(vname)(j,:) = ...
                      oo_recursive_{nobs(j)}.PointForecast.HPDsup.(vname);
            else
                oo_.RecursiveForecast.Mean.(vname)(j,:) =...
                      oo_recursive_{nobs(j)}.forecast.Mean.(vname);
                oo_.RecursiveForecast.HPDinf.(vname)(j,:) =...
                      oo_recursive_{nobs(j)}.forecast.HPDinf.(vname);
                oo_.RecursiveForecast.HPDsup.(vname)(j,:) =...
                      oo_recursive_{nobs(j)}.forecast.HPDsup.(vname);
              end
            x = nobs(1)+nobs(j)-nobs(1)+(1:horizon);

            y = oo_.RecursiveForecast.Mean.(vname)(j,:);
            y1 = oo_.RecursiveForecast.HPDinf.(vname)(j,:);
            y2 = oo_.RecursiveForecast.HPDsup.(vname)(j,:);
            plot(x,y,'-b','linewidth',2)
            plot(x,y1,'--g', ...
                 'linewidth',1.5)
            plot(x,y2,'--g', ...
                 'linewidth',1.5)
            if mh_replic
                y3 = oo_.RecursiveForecast.HPDTotalinf.(vname)(j,:);
                y4 = oo_.RecursiveForecast.HPDTotalsup.(vname)(j,:);
                plot(x,y3,'--r', ...
                     'linewidth',1.5)
                plot(x,y4,'--r','linewidth',1.5)
            end
        end
        box on
        title(vname,'Interpreter','none')
        hold off
        xlim([nobs(1)-offsetx nobs(end)+horizon])
        m = m + 1;
        if mod(i+1,nstar) == 1 || i ==size(var_list,1)
            dyn_saveas(hfig,[M_.fname,filesep,'graphs',filesep M_.fname '_RecursiveForecasts_' int2str(plot_index)],options_);
        end
    end
end
