function oo_ = McMCDiagnostics(options_, estim_params_, M_, oo_)
% function McMCDiagnostics
% Computes convergence tests 
% 
% INPUTS 
%   options_         [structure]
%   estim_params_    [structure]
%   M_               [structure]
%
% OUTPUTS 
%   oo_              [structure] 
%
% SPECIAL REQUIREMENTS
%   none
%
% PARALLEL CONTEXT
% See the comment in posterior_sampler.m funtion.

% Copyright (C) 2005-2016 Dynare Team
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

OutputFolder = CheckPath('Output',M_.dname);
MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;

TeX = options_.TeX;
nblck = options_.mh_nblck;
npar = estim_params_.nvx;
npar = npar + estim_params_.nvn;
npar = npar + estim_params_.ncx;
npar = npar + estim_params_.ncn;
npar = npar + estim_params_.np ;
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);

load_last_mh_history_file(MetropolisFolder, ModelName);
NumberOfMcFilesPerBlock = record.LastFileNumber;

% Check that the declared number of blocks is consistent with informations saved in mh-history files.
if ~isequal(nblck,record.Nblck)
    disp(['Estimation::mcmc::diagnostics: The number of MCMC chains you declared (' num2str(nblck) ') is inconsistent with the information available in the mh-history files (' num2str(record.Nblck) ' chains)!'])
    disp(['                               I reset the number of MCMC chains to ' num2str(record.Nblck) '.'])
    nblck = record.Nblck;
end

% check if all the mh files are available.
issue_an_error_message = 0;
for b = 1:nblck
    nfiles = length(dir([MetropolisFolder ,filesep, ModelName '_mh*_blck' num2str(b) '.mat']));
    if ~isequal(NumberOfMcFilesPerBlock,nfiles)
        issue_an_error_message = 1;
        disp(['Estimation::mcmc::diagnostics: The number of MCMC files in chain ' num2str(b) ' is ' num2str(nfiles) ' while the mh-history files indicate that we should have ' num2str(NumberOfMcFilesPerBlock) ' MCMC files per chain!'])
    end
end
if issue_an_error_message
    error('Estimation::mcmc::diagnostics: I cannot proceed because some MCMC files are missing. Check your MCMC files...')
end

% compute inefficiency factor
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
FirstMhFile = record.KeepedDraws.FirstMhFile;
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);

param_name=[];
param_name_tex=[];
for jj=1:npar
    if options_.TeX
        [par_name_temp,par_name_tex_temp]=get_the_name(jj,options_.TeX,M_,estim_params_,options_);
        param_name = strvcat(param_name,par_name_temp);        
        par_name_tex_temp = strrep(par_name_tex_temp,'$','');
        param_name_tex = strvcat(param_name_tex,par_name_tex_temp);
    else
        [par_name_temp]=get_the_name(jj,options_.TeX,M_,estim_params_,options_);
        param_name = strvcat(param_name,par_name_temp);            
    end
    Draws = GetAllPosteriorDraws(jj,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
    Draws = reshape(Draws,[NumberOfDraws nblck]);
    Nc = min(1000, NumberOfDraws/2);
    for ll=1:nblck
        Ifac(ll,jj) = mcmc_ifac(Draws(:,ll), Nc);
    end
    tmp = num2cell(Ifac(:,jj));
end

my_title='MCMC Inefficiency factors per block';
IFAC_header='Parameter';
IFAC_header_tex='Parameter';
for j=1:nblck,
    IFAC_header = char(IFAC_header,['Block ' int2str(j)]);
    IFAC_header_tex = char(IFAC_header_tex,['Block~' int2str(j)]);
end

lh = size(param_name,2)+2;
dyntable(options_,my_title,IFAC_header,param_name,Ifac',lh,12,3);
skipline()

if options_.TeX
    dyn_latex_table(M_,options_,my_title,'MCMC_inefficiency_factors',IFAC_header_tex,param_name_tex,Ifac',lh,12,3);
end
record.InefficiencyFactorsPerBlock = Ifac;
update_last_mh_history_file(MetropolisFolder, ModelName, record);


PastDraws = sum(record.MhDraws,1);
LastFileNumber = PastDraws(2);
LastLineNumber = record.MhDraws(end,3);
NumberOfDraws  = PastDraws(1);

if NumberOfDraws<=2000
    warning(['estimation:: MCMC convergence diagnostics are not computed because the total number of iterations is less than 2000!'])
    return
end

if nblck == 1 % Brooks and Gelman tests need more than one block
    convergence_diagnostics_geweke=zeros(npar,4+2*length(options_.convergence.geweke.taper_steps));
    if any(options_.convergence.geweke.geweke_interval<0) || any(options_.convergence.geweke.geweke_interval>1) || length(options_.convergence.geweke.geweke_interval)~=2 ...
        || (options_.convergence.geweke.geweke_interval(2)-options_.convergence.geweke.geweke_interval(1)<0)
        fprintf('\nCONVERGENCE DIAGNOSTICS: Invalid option for geweke_interval. Using the default of [0.2 0.5].\n')
        options_.convergence.geweke.geweke_interval=[0.2 0.5];
    end
    first_obs_begin_sample = max(1,ceil(options_.mh_drop*NumberOfDraws));
    last_obs_begin_sample = first_obs_begin_sample+round(options_.convergence.geweke.geweke_interval(1)*NumberOfDraws*(1-options_.mh_drop));
    first_obs_end_sample = first_obs_begin_sample+round(options_.convergence.geweke.geweke_interval(2)*NumberOfDraws*(1-options_.mh_drop));
    param_name=[];
    if options_.TeX
        param_name_tex=[];    
    end
    for jj=1:npar        
        if options_.TeX
            [param_name_temp, param_name_tex_temp]= get_the_name(jj,options_.TeX,M_,estim_params_,options_);
            param_name_tex = strvcat(param_name_tex,strrep(param_name_tex_temp,'$',''));
            param_name = strvcat(param_name,param_name_temp);
        else
            param_name_temp = get_the_name(jj,options_.TeX,M_,estim_params_,options_);
            param_name = strvcat(param_name,param_name_temp);
        end
    end
    fprintf('\nGeweke (1992) Convergence Tests, based on means of draws %d to %d vs %d to %d.\n',first_obs_begin_sample,last_obs_begin_sample,first_obs_end_sample,NumberOfDraws);
    fprintf('p-values are for Chi2-test for equality of means.\n');    
    Geweke_header=char('Parameter', 'Post. Mean', 'Post. Std', 'p-val No Taper');
    for ii=1:length(options_.convergence.geweke.taper_steps)
        Geweke_header=char(Geweke_header,['p-val ' num2str(options_.convergence.geweke.taper_steps(ii)),'% Taper']);
    end
    datamat=NaN(npar,3+length(options_.convergence.geweke.taper_steps));
    for jj=1:npar
        startline=0;
        for n = 1:NumberOfMcFilesPerBlock
            load([MetropolisFolder '/' ModelName '_mh',int2str(n),'_blck1.mat'],'x2');
            nx2 = size(x2,1);
            param_draws(startline+(1:nx2),1) = x2(:,jj);
            startline = startline + nx2;
        end
        [results_vec, results_struct] = geweke_moments(param_draws,options_);
        convergence_diagnostics_geweke(jj,:)=results_vec;
    
        param_draws1 = param_draws(first_obs_begin_sample:last_obs_begin_sample,:);
        param_draws2 = param_draws(first_obs_end_sample:end,:);
        [results_vec1] = geweke_moments(param_draws1,options_);
        [results_vec2] = geweke_moments(param_draws2,options_);
        
        results_struct = geweke_chi2_test(results_vec1,results_vec2,results_struct,options_);
        eval(['oo_.convergence.geweke.',param_name(jj,:),'=results_struct;'])
        datamat(jj,:)=[results_struct.posteriormean,results_struct.posteriorstd,results_struct.prob_chi2_test];
    end
    lh = size(param_name,2)+2;
    dyntable(options_,'',Geweke_header,param_name,datamat,lh,12,3);
    if options_.TeX
        Geweke_tex_header=char('Parameter', 'Mean', 'Std', 'No\ Taper');
        additional_header={[' & \multicolumn{2}{c}{Posterior} & \multicolumn{',num2str(1+length(options_.convergence.geweke.taper_steps)),'}{c}{p-values} \\'],
            ['\cmidrule(r{.75em}){2-3} \cmidrule(r{.75em}){4-',num2str(4+length(options_.convergence.geweke.taper_steps)),'}']};
        for ii=1:length(options_.convergence.geweke.taper_steps)
            Geweke_tex_header=char(Geweke_tex_header,[num2str(options_.convergence.geweke.taper_steps(ii)),'\%%\ Taper']);
        end
        headers = char(Geweke_tex_header);
        lh = size(param_name_tex,2)+2;
        my_title=sprintf('Geweke (1992) Convergence Tests, based on means of draws %d to %d vs %d to %d. p-values are for $\\\\chi^2$-test for equality of means.',first_obs_begin_sample,last_obs_begin_sample,first_obs_end_sample,NumberOfDraws);
        dyn_latex_table(M_,options_,my_title,'geweke',headers,param_name_tex,datamat,lh,12,4,additional_header);    
    end
    skipline(2);
    return;
end

Origin = 1000;
StepSize = ceil((NumberOfDraws-Origin)/100);% So that the computational time does not 
ALPHA = 0.2;                                % increase too much with the number of simulations. 
time = 1:NumberOfDraws;
xx = Origin:StepSize:NumberOfDraws;
NumberOfLines = length(xx);
tmp = zeros(NumberOfDraws*nblck,3);
UDIAG = zeros(NumberOfLines,6,npar);

if NumberOfDraws < Origin
    disp('Estimation::mcmc::diagnostics: The number of simulations is too small to compute the MCMC convergence diagnostics.')
    return
end

if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([OutputFolder '/' ModelName '_UnivariateDiagnostics.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by McmcDiagnostics.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

disp('Estimation::mcmc::diagnostics: Univariate convergence diagnostic, Brooks and Gelman (1998):')

% The mandatory variables for local/remote parallel
% computing are stored in localVars struct.

localVars.MetropolisFolder = MetropolisFolder;
localVars.nblck = nblck;
localVars.NumberOfMcFilesPerBlock = NumberOfMcFilesPerBlock;
localVars.Origin = Origin;
localVars.StepSize = StepSize;
localVars.mh_drop = options_.mh_drop;
localVars.NumberOfDraws = NumberOfDraws;
localVars.NumberOfLines = NumberOfLines;
localVars.time = time;
localVars.M_ = M_;


% Like sequential execution!
if isnumeric(options_.parallel),
    fout = McMCDiagnostics_core(localVars,1,npar,0);
    UDIAG = fout.UDIAG;
    clear fout
    % Parallel execution!
else
    ModelName = ModelName;
    if ~isempty(M_.bvar)
        ModelName = [ModelName '_bvar'];
    end
    NamFileInput={[M_.dname '/metropolis/'],[ModelName '_mh*_blck*.mat']};
    
    [fout, nBlockPerCPU, totCPU] = masterParallel(options_.parallel, 1, npar,NamFileInput,'McMCDiagnostics_core', localVars, [], options_.parallel_info);
    UDIAG = fout(1).UDIAG;
    for j=2:totCPU,
        UDIAG = cat(3,UDIAG ,fout(j).UDIAG);
    end
end

UDIAG(:,[2 4 6],:) = UDIAG(:,[2 4 6],:)/nblck;
skipline()
clear pmet temp moyenne CSUP CINF csup cinf n linea iter tmp;    
pages = floor(npar/3);
k = 0;  
for i = 1:pages
    h=dyn_figure(options_,'Name','MCMC univariate convergence diagnostic (Brooks and Gelman,1998)');
    boxplot = 1;
    for j = 1:3 % Loop over parameters
        k = k+1;
        [nam,namtex] = get_the_name(k,TeX,M_,estim_params_,options_);
        for crit = 1:3% Loop over criteria
            if crit == 1
                plt1 = UDIAG(:,1,k);
                plt2 = UDIAG(:,2,k);
                namnam  = [nam , ' (Interval)']; 
            elseif crit == 2
                plt1 = UDIAG(:,3,k);
                plt2 = UDIAG(:,4,k);
                namnam  = [nam , ' (m2)'];
            elseif crit == 3    
                plt1 = UDIAG(:,5,k);
                plt2 = UDIAG(:,6,k);
                namnam  = [nam , ' (m3)'];
            end
            if TeX
                if j==1
                    NAMES = deblank(namnam);
                    TEXNAMES = deblank(namtex);
                else
                    NAMES = char(NAMES,deblank(namnam));
                    TEXNAMES = char(TEXNAMES,deblank(namtex));
                end
            end
            subplot(3,3,boxplot);
            plot(xx,plt1,'-b');     % Pooled
            hold on;
            plot(xx,plt2,'-r');     % Within (mean)
            hold off;
            xlim([xx(1) xx(NumberOfLines)])
            title(namnam,'Interpreter','none')
            boxplot = boxplot + 1;
        end
    end
    dyn_saveas(h,[OutputFolder '/' ModelName '_udiag' int2str(i)],options_);
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:size(NAMES,1)
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_udiag%s}\n',options_.figures.textwidth*min((boxplot-1)/3,1),[OutputFolder '/' ModelName],int2str(i));
        fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
        fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
        fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
        fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(i));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
    end
end
reste = npar-k;
if reste
    if reste == 1
        nr = 3;
        nc = 1;
    elseif reste == 2;
        nr = 2;
        nc = 3;
    end
    h = dyn_figure(options_,'Name','MCMC univariate convergence diagnostic (Brooks and Gelman, 1998)');
    boxplot = 1;
    for j = 1:reste
        k = k+1;
        [nam,namtex] = get_the_name(k,TeX,M_,estim_params_,options_);
        for crit = 1:3
            if crit == 1
                plt1 = UDIAG(:,1,k);
                plt2 = UDIAG(:,2,k);
                namnam  = [nam , ' (Interval)']; 
            elseif crit == 2
                plt1 = UDIAG(:,3,k);
                plt2 = UDIAG(:,4,k);
                namnam  = [nam , ' (m2)'];
            elseif crit == 3    
                plt1 = UDIAG(:,5,k);
                plt2 = UDIAG(:,6,k);
                namnam  = [nam , ' (m3)'];
            end
            if TeX
                if j==1
                    NAMES = deblank(namnam);
                    TEXNAMES = deblank(namtex);
                else
                    NAMES = char(NAMES,deblank(namnam));
                    TEXNAMES = char(TEXNAMES,deblank(namtex));
                end
            end
            subplot(nr,nc,boxplot);
            plot(xx,plt1,'-b');                                 % Pooled
            hold on;
            plot(xx,plt2,'-r');                                 % Within (mean)
            hold off;
            xlim([xx(1) xx(NumberOfLines)]);
            title(namnam,'Interpreter','none');
            boxplot = boxplot + 1;
        end
    end
    dyn_saveas(h,[ OutputFolder '/' ModelName '_udiag' int2str(pages+1)],options_);
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:size(NAMES,1);
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_udiag%s}\n',options_.figures.textwidth*min((boxplot-1)/nc,1),[OutputFolder '/' ModelName],int2str(pages+1));
        if reste == 2
            fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
            fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
            fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
        elseif reste == 1
            fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
            fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
            fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
        end
        fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(pages+1));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'% End Of TeX file.');
        fclose(fidTeX);
    end
end % if reste > 0
clear UDIAG;
%
% Multivariate diagnostic.
%
if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([OutputFolder '/' ModelName '_MultivariateDiagnostics.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by McmcDiagnostics.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end
tmp = zeros(NumberOfDraws*nblck,3);
MDIAG = zeros(NumberOfLines,6);
for b = 1:nblck
    startline = 0;
    for n = 1:NumberOfMcFilesPerBlock
        load([MetropolisFolder '/' ModelName '_mh',int2str(n),'_blck' int2str(b) '.mat'],'logpo2');
        nlogpo2 = size(logpo2,1);
        tmp((b-1)*NumberOfDraws+startline+(1:nlogpo2),1) = logpo2;
        startline = startline+nlogpo2;
    end
end
clear logpo2;
tmp(:,2) = kron(transpose(1:nblck),ones(NumberOfDraws,1));
tmp(:,3) = kron(ones(nblck,1),time'); 
tmp = sortrows(tmp,1);
ligne   = 0;
for iter  = Origin:StepSize:NumberOfDraws
    ligne = ligne+1;
    linea = ceil(options_.mh_drop*iter);
    n     = iter-linea+1;
    cinf  = round(n*ALPHA/2);
    csup  = round(n*(1-ALPHA/2));
    CINF  = round(nblck*n*ALPHA/2);
    CSUP  = round(nblck*n*(1-ALPHA/2));
    temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
    MDIAG(ligne,1) = temp(CSUP,1)-temp(CINF,1);
    moyenne = mean(temp(:,1));%% Pooled mean.
    MDIAG(ligne,3) = sum((temp(:,1)-moyenne).^2)/(nblck*n-1);
    MDIAG(ligne,5) = sum(abs(temp(:,1)-moyenne).^3)/(nblck*n-1);
    for i=1:nblck
        pmet = temp(find(temp(:,2)==i));
        MDIAG(ligne,2) = MDIAG(ligne,2) + pmet(csup,1)-pmet(cinf,1);
        moyenne = mean(pmet,1); %% Within mean. 
        MDIAG(ligne,4) = MDIAG(ligne,4) + sum((pmet(:,1)-moyenne).^2)/(n-1);
        MDIAG(ligne,6) = MDIAG(ligne,6) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
    end
end
MDIAG(:,[2 4 6],:) = MDIAG(:,[2 4 6],:)/nblck;  

h = dyn_figure(options_,'Name','Multivariate convergence diagnostic');
boxplot = 1;
for crit = 1:3
    if crit == 1
        plt1 = MDIAG(:,1);
        plt2 = MDIAG(:,2);
        namnam  = 'Interval'; 
    elseif crit == 2
        plt1 = MDIAG(:,3);
        plt2 = MDIAG(:,4);
        namnam  = 'm2';
    elseif crit == 3    
        plt1 = MDIAG(:,5);
        plt2 = MDIAG(:,6);
        namnam  = 'm3';
    end
    if TeX
        if crit == 1
            NAMES = deblank(namnam);
        else
            NAMES = char(NAMES,deblank(namnam));
        end
    end
    subplot(3,1,boxplot);
    plot(xx,plt1,'-b');  % Pooled
    hold on
    plot(xx,plt2,'-r');  % Within (mean)
    hold off
    xlim([xx(1) xx(NumberOfLines)])
    title(namnam,'Interpreter','none');
    boxplot = boxplot + 1;
end
dyn_saveas(h,[ OutputFolder '/' ModelName '_mdiag'],options_);

if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    for jj = 1:3
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),' ');
    end    
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_mdiag}\n',[OutputFolder '/' ModelName]);
    fprintf(fidTeX,'\\caption{Multivariate convergence diagnostics for the Metropolis-Hastings.\n');
    fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
    fprintf(fidTeX,'the eighty percent interval, the second and third moments. The different \n');
    fprintf(fidTeX,'parameters are aggregated using the posterior kernel.}');
    fprintf(fidTeX,'\\label{Fig:MultivariateDiagnostics}\n');
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,'\n');
    fprintf(fidTeX,'%% End Of TeX file.');
    fclose(fidTeX);
end

