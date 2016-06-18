function pm3(n1,n2,ifil,B,tit1,tit2,tit3,tit_tex,names1,names2,name3,DirectoryName,var_type)
% pm3(n1,n2,ifil,B,tit1,tit2,tit3,tit_tex,names1,names2,name3,DirectoryName,var_type)
% Computes, stores and plots the posterior moment statistics.
%
% Inputs:
%  n1           [scalar] size of first dimension of moment matrix
%  n2           [scalar] size of second dimension of moment matrix
%  ifil         [scalar] number of moment files to load
%  B            [scalar] number of subdraws
%  tit1         [string] Figure title 
%  tit2         [string] not used
%  tit3         [string] Save name for figure
%  tit_tex      [cell array] TeX-Names for Variables
%  names1       [cell array] Names of all variables in the moment matrix from
%                       which names2 is selected
%  names2       [cell array] Names of variables subset selected for moments
%  names3       [string] Name of the field in oo_ structure to be set
%  DirectoryName [string] Name of the directory in which to save and from
%                       where to read
%  var_type     [string] suffix of the filename from which to load moment
%                   matrix

% PARALLEL CONTEXT
% See also the comment in posterior_sampler.m funtion.


% Copyright (C) 2007-2016 Dynare Team
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

global options_ M_ oo_

nn = 3;
MaxNumberOfPlotsPerFigure = nn^2; % must be square
varlist = names2;
if isempty(varlist)
    varlist = names1;
    SelecVariables = (1:M_.endo_nbr)';
    nvar = M_.endo_nbr;
else
    nvar = size(varlist,1);
    SelecVariables = [];
    for i=1:nvar
        if ~isempty(strmatch(varlist(i,:),names1,'exact'))
            SelecVariables = [SelecVariables;strmatch(varlist(i,:),names1,'exact')];
        end
    end
end
if options_.TeX
    if isempty(tit_tex),
        tit_tex=names1;
    end
        
    varlist_TeX = [];
    for i=1:nvar
        if i==1
            varlist_TeX = tit_tex(SelecVariables(i),:);
        else
            varlist_TeX = char(varlist_TeX,tit_tex(SelecVariables(i),:));
        end
    end
end
Mean = zeros(n2,nvar);
Median = zeros(n2,nvar);
Var = zeros(n2,nvar);
Distrib = zeros(9,n2,nvar);
HPD = zeros(2,n2,nvar);
if options_.estimation.moments_posterior_density.indicator
    Density = zeros(options_.estimation.moments_posterior_density.gridpoints,2,n2,nvar);
end
fprintf(['Estimation::mcmc: ' tit1 '\n']);
k = 0;
filter_step_ahead_indicator=0;
filter_covar_indicator=0;

for file = 1:ifil
    load([DirectoryName '/' M_.fname var_type int2str(file)]);
    if strcmp(var_type,'_filter_step_ahead')
        if file==1 %on first run, initialize variable for storing filter_step_ahead
            stock1_filter_step_ahead=NaN(n1,n2,B,length(options_.filter_step_ahead)); 
            stock1 = zeros(n1,n2,B);
        end
        filter_step_ahead_indicator=1;
        stock_filter_step_ahead=zeros(n1,n2,size(stock,4),length(options_.filter_step_ahead));
        for ii=1:length(options_.filter_step_ahead)
            K_step_ahead=options_.filter_step_ahead(ii);
            stock_filter_step_ahead(:,:,:,ii)=stock(ii,:,1+K_step_ahead:n2+K_step_ahead,:);
        end
        stock = squeeze(stock(1,:,1+1:1+n2,:)); %1 step ahead starts at entry 2
        k = k(end)+(1:size(stock,3));
        stock1(:,:,k) = stock;
        stock1_filter_step_ahead(:,:,k,:) = stock_filter_step_ahead;
    elseif strcmp(var_type,'_filter_covar')
        if file==1 %on first run, initialize variable for storing filter_step_ahead
            stock1_filter_covar=NaN(n1,n2,size(stock,3),B);
        end
        filter_covar_indicator=1;
        k = k(end)+(1:size(stock,4));
        stock1_filter_covar(:,:,:,k) = stock;
    elseif strcmp(var_type,'_trend_coeff')
        if file==1 %on first run, initialize variable for storing filter_step_ahead
            stock1_filter_step_ahead=NaN(n1,n2,B,length(options_.filter_step_ahead));
            stock1 = zeros(n1,B);
        end
        k = k(end)+(1:size(stock,2));
        stock1(:,k) = stock;
    else
        if file==1 %on first run, initialize variable for storing filter_step_ahead
            stock1 = zeros(n1,n2,B);
        end
        k = k(end)+(1:size(stock,3));
        stock1(:,:,k) = stock;
    end
end
clear stock
if filter_step_ahead_indicator
    clear stock_filter_step_ahead
    filter_steps=length(options_.filter_step_ahead);
    Mean_filter_step_ahead = zeros(filter_steps,nvar,n2);
    Median_filter_step_ahead = zeros(filter_steps,nvar,n2);
    Var_filter_step_ahead = zeros(filter_steps,nvar,n2);
    Distrib_filter_step_ahead = zeros(9,filter_steps,nvar,n2);
    HPD_filter_step_ahead = zeros(2,filter_steps,nvar,n2);
    if options_.estimation.moments_posterior_density.indicator
        Density_filter_step_ahead = zeros(options_.estimation.moments_posterior_density.gridpoints,2,filter_steps,nvar,n2);
    end
end
if filter_covar_indicator
    draw_dimension=4;
    oo_.FilterCovariance.Mean = squeeze(mean(stock1_filter_covar,draw_dimension));
    oo_.FilterCovariance.Median = squeeze(median(stock1_filter_covar,draw_dimension));
    oo_.FilterCovariance.var = squeeze(var(stock1_filter_covar,0,draw_dimension));
    if size(stock1_filter_covar,draw_dimension)>2
        hpd_interval = quantile(stock1_filter_covar,[(1-options_.mh_conf_sig)/2 (1-options_.mh_conf_sig)/2+options_.mh_conf_sig],draw_dimension);
    else
        size_matrix=size(stock1_filter_covar);
        hpd_interval=NaN([size_matrix(1:3),2]);
    end
    if size(stock1_filter_covar,draw_dimension)>9
        post_deciles =quantile(stock1_filter_covar,[0.1:0.1:0.9],draw_dimension);
    else
        size_matrix=size(stock1_filter_covar);
        post_deciles=NaN([size_matrix(1:3),9]);
    end
    oo_.FilterCovariance.post_deciles=post_deciles;
    oo_.FilterCovariance.HPDinf=squeeze(hpd_interval(:,:,:,1));
    oo_.FilterCovariance.HPDsup=squeeze(hpd_interval(:,:,:,2));
    fprintf(['Estimation::mcmc: ' tit1 ', done!\n']);
    return
end

if strcmp(var_type,'_trend_coeff') %two dimensional arrays
    for i = 1:nvar
        if options_.estimation.moments_posterior_density.indicator
            [Mean(1,i),Median(1,i),Var(1,i),HPD(:,1,i),Distrib(:,1,i),Density(:,:,1,i)] = ...
                posterior_moments(squeeze(stock1(SelecVariables(i),:)),1,options_.mh_conf_sig,options_.estimation.moments_posterior_density);
        else
            [Mean(1,i),Median(1,i),Var(1,i),HPD(:,1,i),Distrib(:,1,i)] = ...
                posterior_moments(squeeze(stock1(SelecVariables(i),:)),0,options_.mh_conf_sig);
        end
    end
else %three dimensional arrays
    for i = 1:nvar
        for j = 1:n2
            if options_.estimation.moments_posterior_density.indicator
                [Mean(j,i),Median(j,i),Var(j,i),HPD(:,j,i),Distrib(:,j,i),Density(:,:,j,i)] = ...
                    posterior_moments(squeeze(stock1(SelecVariables(i),j,:)),1,options_.mh_conf_sig,options_.estimation.moments_posterior_density);
            else
                [Mean(j,i),Median(j,i),Var(j,i),HPD(:,j,i),Distrib(:,j,i)] = ...
                    posterior_moments(squeeze(stock1(SelecVariables(i),j,:)),0,options_.mh_conf_sig);
            end
            if filter_step_ahead_indicator
                if options_.estimation.moments_posterior_density.indicator
                    for K_step = 1:length(options_.filter_step_ahead)
                        [Mean_filter_step_ahead(K_step,i,j),Median_filter_step_ahead(K_step,i,j),Var_filter_step_ahead(K_step,i,j),HPD_filter_step_ahead(:,K_step,i,j),Distrib_filter_step_ahead(:,K_step,i,j),Density_filter_step_ahead(:,:,K_step,i,j) ] = ...
                            posterior_moments(squeeze(stock1_filter_step_ahead(SelecVariables(i),j,:,K_step)),1,options_.mh_conf_sig,options_.estimation.moments_posterior_density);
                    end
                else
                    for K_step = 1:length(options_.filter_step_ahead)
                        [Mean_filter_step_ahead(K_step,i,j),Median_filter_step_ahead(K_step,i,j),Var_filter_step_ahead(K_step,i,j),HPD_filter_step_ahead(:,K_step,i,j),Distrib_filter_step_ahead(:,K_step,i,j)] = ...
                            posterior_moments(squeeze(stock1_filter_step_ahead(SelecVariables(i),j,:,K_step)),0,options_.mh_conf_sig);
                    end
                end
            end
        end
    end
end

clear stock1
if filter_step_ahead_indicator %write matrices corresponding to ML
    clear stock1_filter_step_ahead
    FilteredVariablesKStepAhead=zeros(length(options_.filter_step_ahead),nvar,n2+max(options_.filter_step_ahead));
    FilteredVariablesKStepAheadVariances=zeros(length(options_.filter_step_ahead),nvar,n2+max(options_.filter_step_ahead));
    for K_step = 1:length(options_.filter_step_ahead)
        FilteredVariablesKStepAhead(K_step,:,1+options_.filter_step_ahead(K_step):n2+options_.filter_step_ahead(K_step))=Mean_filter_step_ahead(K_step,:,:);
        FilteredVariablesKStepAheadVariances(K_step,:,1+options_.filter_step_ahead(K_step):n2+options_.filter_step_ahead(K_step))=Mean_filter_step_ahead(K_step,:,:);
    end
    oo_.FilteredVariablesKStepAhead=FilteredVariablesKStepAhead;
    oo_.FilteredVariablesKStepAheadVariances=FilteredVariablesKStepAheadVariances;
end

if strcmp(var_type,'_trend_coeff') || strcmp(var_type,'_smoothed_trend') || strcmp(var_type,'_smoothed_trend')
    for i = 1:nvar
        name = deblank(names1(SelecVariables(i),:));
        oo_.Smoother.(name3).Mean.(name) = Mean(:,i);
        oo_.Smoother.(name3).Median.(name) = Median(:,i);
        oo_.Smoother.(name3).Var.(name) = Var(:,i);
        oo_.Smoother.(name3).deciles.(name) = Distrib(:,:,i);
        oo_.Smoother.(name3).HPDinf.(name) = HPD(1,:,i)';
        oo_.Smoother.(name3).HPDsup.(name) = HPD(2,:,i)';
        if options_.estimation.moments_posterior_density.indicator
            oo_.Smoother.(name3).density.(name) = Density(:,:,:,i);
        end
    end
else
    for i = 1:nvar
        name = deblank(names1(SelecVariables(i),:));
        oo_.(name3).Mean.(name) = Mean(:,i);
        oo_.(name3).Median.(name) = Median(:,i);
        oo_.(name3).Var.(name) = Var(:,i);
        oo_.(name3).deciles.(name) = Distrib(:,:,i);
        oo_.(name3).HPDinf.(name) = HPD(1,:,i)';
        oo_.(name3).HPDsup.(name) = HPD(2,:,i)';
        if options_.estimation.moments_posterior_density.indicator
            oo_.(name3).density.(name) = Density(:,:,:,i);
        end
        if filter_step_ahead_indicator
            for K_step = 1:length(options_.filter_step_ahead)
                name4=['Filtered_Variables_',num2str(K_step),'_step_ahead'];
                oo_.(name4).Mean.(name) = squeeze(Mean_filter_step_ahead(K_step,i,:));
                oo_.(name4).Median.(name) = squeeze(Median_filter_step_ahead(K_step,i,:));
                oo_.(name4).Var.(name) = squeeze(Var_filter_step_ahead(K_step,i,:));
                oo_.(name4).deciles.(name) = squeeze(Distrib_filter_step_ahead(:,K_step,i,:));
                oo_.(name4).HPDinf.(name) = squeeze(HPD_filter_step_ahead(1,K_step,i,:));
                oo_.(name4).HPDsup.(name) = squeeze(HPD_filter_step_ahead(2,K_step,i,:));
                if options_.estimation.moments_posterior_density.indicator
                    oo_.(name4).density.(name) = squeeze(Density_filter_step_ahead(:,:,K_step,i,:));
                end
            end
        end
    end
end

if strcmp(var_type,'_trend_coeff') || max(max(abs(Mean(:,:))))<=10^(-6) || all(all(isnan(Mean)))
    fprintf(['Estimation::mcmc: ' tit1 ', done!\n']);
    return %not do plots
end
%%
%% 	Finally I build the plots.
%%

% Block of code executed in parallel, with the exception of file
% .tex generation always run sequentially. This portion of code is execute in parallel by
% pm3_core1.m function.

% %%%%%%%%%   PARALLEL BLOCK % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%% The file .TeX! are not saved in parallel.



% Store the variable mandatory for local/remote parallel computing.

localVars=[];

localVars.tit1=tit1;
localVars.nn=nn;
localVars.n2=n2;
localVars.Distrib=Distrib;
localVars.varlist=varlist;
localVars.MaxNumberOfPlotsPerFigure=MaxNumberOfPlotsPerFigure;
localVars.name3=name3;
localVars.tit3=tit3;
localVars.Mean=Mean;
% Like sequential execution!
nvar0=nvar;

if ~isoctave
    % Commenting for testing!
    if isnumeric(options_.parallel) || ceil(size(varlist,1)/MaxNumberOfPlotsPerFigure)<4,
        fout = pm3_core(localVars,1,nvar,0);
        
        % Parallel execution!
    else
        isRemoteOctave = 0;
        for indPC=1:length(options_.parallel),
            isRemoteOctave = isRemoteOctave + (findstr(options_.parallel(indPC).MatlabOctavePath, 'octave'));
        end
        if isRemoteOctave
            fout = pm3_core(localVars,1,nvar,0);
        else
            globalVars = struct('M_',M_, ...
                'options_', options_, ...
                'oo_', oo_);
            [fout, nvar0, totCPU] = masterParallel(options_.parallel, 1, nvar, [],'pm3_core', localVars,globalVars, options_.parallel_info);
        end
    end
else
    % For the time being in Octave enviroment the pm3.m is executed only in
    % serial modality, to avoid problem with the plots.
    
    fout = pm3_core(localVars,1,nvar,0);
end

subplotnum = 0;

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([M_.dname '/Output/' M_.fname '_' name3 '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by Dynare.\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
    nvar0=cumsum(nvar0);

    i=0;    
    for j=1:length(nvar0),
    
    NAMES = [];
    TEXNAMES = [];
    nvar=nvar0(j);
    while i<nvar,
        i=i+1;
        if max(abs(Mean(:,i))) > 10^(-6)
            subplotnum = subplotnum+1;
            name = deblank(varlist(i,:));
            texname = deblank(varlist_TeX(i,:));
            if subplotnum==1
                NAMES = name;
                TEXNAMES = ['$' texname '$'];
            else
                NAMES = char(NAMES,name);
                TEXNAMES = char(TEXNAMES,['$' texname '$']);
            end
        end
        if subplotnum == MaxNumberOfPlotsPerFigure || i == nvar
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:size(TEXNAMES,1)
                fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
            end
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.8\\textwidth]{%s/Output/%s_' name3 '_%s}\n'],M_.dname,M_.fname,deblank(tit3(i,:)));
            fprintf(fidTeX,'\\label{Fig:%s:%s}\n',name3,deblank(tit3(i,:)));
            fprintf(fidTeX,'\\caption{%s}\n',tit1);
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            subplotnum = 0;
            NAMES = [];
            TEXNAMES = [];
        end
    end
    end
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
end

fprintf(['Estimation::mcmc: ' tit1 ', done!\n']);