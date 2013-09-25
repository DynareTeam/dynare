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
% See the comment in random_walk_metropolis_hastings.m funtion.

% Copyright (C) 2005-2013 Dynare Team
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

DirectoryName = CheckPath('Output',M_.dname);
MhDirectoryName = CheckPath('metropolis',M_.dname);

TeX = options_.TeX;
nblck = options_.mh_nblck;
npar = estim_params_.nvx;
npar = npar + estim_params_.nvn;
npar = npar + estim_params_.ncx;
npar = npar + estim_params_.ncn;
npar = npar + estim_params_.np ;
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);

load([MhDirectoryName '/'  M_.fname '_mh_history.mat'])

NumberOfMcFilesPerBlock = size(dir([MhDirectoryName ,filesep, M_.fname '_mh*_blck1.mat']),1);

% check if all previous files are there for block 1
check_presence_consecutive_MC_files(MhDirectoryName,M_.fname,1)

if nblck == 1 % Brooks and Gelman tests need more than one block 
    convergence_diagnostics_geweke=zeros(npar,4+2*length(options_.convergence.geweke.taper_steps));
    if any(options_.convergence.geweke.geweke_interval<0) || any(options_.convergence.geweke.geweke_interval>1) || length(any(options_.convergence.geweke.geweke_interval<0))~=2 ...
        || (options_.convergence.geweke.geweke_interval(2)-options_.convergence.geweke.geweke_interval(1)<0)
        fprintf('\nCONVERGENCE DIAGNOSTICS: Invalid option for geweke_interval. Using the default of [0.2 0.5].\n')
        options_.convergence.geweke.geweke_interval=[0.2 0.5];
    end
    first_obs_begin_sample = max(1,ceil(options_.mh_drop*options_.mh_replic));
    last_obs_begin_sample = first_obs_begin_sample+round(options_.convergence.geweke.geweke_interval(1)*options_.mh_replic*options_.mh_drop);
    first_obs_end_sample = first_obs_begin_sample+round(options_.convergence.geweke.geweke_interval(2)*options_.mh_replic*options_.mh_drop);
    param_name=[];
    for jj=1:npar
        param_name = strvcat(param_name,get_the_name(jj,options_.TeX,M_,estim_params_,options_));
    end
    fprintf('\nGeweke (1992) Convergence Tests, based on means of draws %d to %d vs %d to %d.\n',first_obs_begin_sample,last_obs_begin_sample,first_obs_end_sample,options_.mh_replic);
    fprintf('p-values are for Chi2-test for equality of means.\n');    
    Geweke_header={'Parameter', 'Post. Mean', 'Post. Std', 'p-val No Taper'};
    print_string=['%',num2str(size(param_name,2)+3),'s \t %12.3f \t %12.3f \t %12.3f'];
    print_string_header=['%',num2str(size(param_name,2)+3),'s \t %12s \t %12s \t %12s'];    
    for ii=1:length(options_.convergence.geweke.taper_steps)
        Geweke_header=[Geweke_header, ['p-val ' num2str(options_.convergence.geweke.taper_steps(ii)),'% Taper']];
        print_string=[print_string,'\t %12.3f'];
        print_string_header=[print_string_header,'\t %12s'];
    end
    print_string=[print_string,'\n'];
    print_string_header=[print_string_header,'\n'];
    fprintf(print_string_header,Geweke_header{1,:});
    for jj=1:npar
        startline=0;
        for n = 1:NumberOfMcFilesPerBlock
            load([MhDirectoryName '/' M_.fname '_mh',int2str(n),'_blck1.mat'],'x2');
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
        fprintf(print_string,param_name(jj,:),results_struct.posteriormean,results_struct.posteriorstd,results_struct.prob_chi2_test)
    end
    skipline(2);
    return;
end

for blck = 2:nblck
    tmp = size(dir([MhDirectoryName ,filesep, M_.fname '_mh*_blck' int2str(blck) '.mat']),1);
    if tmp~=NumberOfMcFilesPerBlock
        disp(['McMCDiagnostics:: The number of mh files in chain ' int2str(blck) ' is ' int2str(tmp) ' while'])
        disp(['                  the number of mh files in chain 1 is ' int2str(mcfiles) '!'])
        error('The number of mh files has to be constant across chains!')
    end
    % check if all previous files are there for block blck
    check_presence_consecutive_MC_files(MhDirectoryName,M_.fname,blck)
end

PastDraws = sum(record.MhDraws,1);
LastFileNumber = PastDraws(2);
LastLineNumber = record.MhDraws(end,3);
NumberOfDraws  = PastDraws(1);

Origin = 1000;
StepSize = ceil((NumberOfDraws-Origin)/100);% So that the computational time does not 
ALPHA = 0.2;                                % increase too much with the number of simulations. 
time = 1:NumberOfDraws;
xx = Origin:StepSize:NumberOfDraws;
NumberOfLines = length(xx);
tmp = zeros(NumberOfDraws*nblck,3);
UDIAG = zeros(NumberOfLines,6,npar);

if NumberOfDraws < Origin
    disp('MCMC Diagnostics :: The number of simulations is to small to compute the MCMC convergence diagnostics.')
    return
end

if TeX
    fidTeX = fopen([DirectoryName '/' M_.fname '_UnivariateDiagnostics.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by McmcDiagnostics.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

disp('MCMC Diagnostics: Univariate convergence diagnostic, Brooks and Gelman (1998):')

% The mandatory variables for local/remote parallel
% computing are stored in localVars struct.

localVars.MhDirectoryName = MhDirectoryName;
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
    ModelName = M_.fname;
    if ~isempty(M_.bvar)
        ModelName = [M_.fname '_bvar'];
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
    dyn_saveas(h,[DirectoryName '/' M_.fname '_udiag' int2str(i)],options_);
    if TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:size(NAMES,1)
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_udiag%s}\n',[DirectoryName '/' M_.fname],int2str(i));
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
    dyn_saveas(h,[ DirectoryName '/' M_.fname '_udiag' int2str(pages+1)],options_);
    if TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:size(NAMES,1);
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_udiag%s}\n',[DirectoryName '/' M_.fname],int2str(pages+1));
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
%%
%% Multivariate diagnostic.
%%
if TeX
    fidTeX = fopen([DirectoryName '/' M_.fname '_MultivariateDiagnostics.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by McmcDiagnostics.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end
tmp = zeros(NumberOfDraws*nblck,3);
MDIAG = zeros(NumberOfLines,6);
for b = 1:nblck
    startline = 0;
    for n = 1:NumberOfMcFilesPerBlock
        %load([MhDirectoryName '/' mcfiles(n,1,b).name],'logpo2');
        load([MhDirectoryName '/' M_.fname '_mh',int2str(n),'_blck' int2str(b) '.mat'],'logpo2');
        nlogpo2 = size(logpo2,1);
        tmp((b-1)*NumberOfDraws+startline+(1:nlogpo2),1) = logpo2;
        startline = startline+nlogpo2;
    end
% $$$   %load([MhDirectoryName '/' mcfiles(NumberOfMcFilesPerBlock,1,b).name],'logpo2');
% $$$   load([MhDirectoryName '/' M_.fname '_mh',int2str(NumberOfMcFilesPerBlock),'_blck' int2str(b) '.mat'],'logpo2');
% $$$   tmp((b-1)*NumberOfDraws+startline+1:(b-1)*NumberOfDraws+ MAX_nruns*(LastFileNumber-1)+LastLineNumber,1) = logpo2;
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
dyn_saveas(h,[ DirectoryName '/' M_.fname '_mdiag'],options_);

if TeX
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    for jj = 1:3
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),' ');
    end    
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_mdiag}\n',[DirectoryName '/' M_.fname]);
    fprintf(fidTeX,'\\caption{Multivariate convergence diagnostics for the Metropolis-Hastings.\n');
    fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
    fprintf(fidTeX,'the eighty percent interval, the second and third moments. The different \n');
    fprintf(fidTeX,'parameters are aggregated using the posterior kernel.}');
    fprintf(fidTeX,'\\label{Fig:MultivariateDiagnostics}\n');
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,'\n');
    fprintf(fidTeX,'% End Of TeX file.');
    fclose(fidTeX);
end

function check_presence_consecutive_MC_files(MhDirectoryName,fname,blck)
% check if all previous files are there
files=ls([MhDirectoryName ,filesep, fname '_mh*_blck' int2str(blck) '.mat']); %list files
right_string=files(:,length([fname '_mh'])+1:end); %cut off left part of filename
k = cell2mat(strfind(cellstr(right_string),['_blck' int2str(blck) '.mat'])); %find index of position after number
file_numbers=str2num(right_string(:,1:k-1)); %get file number
if ~isempty(file_numbers) && ...
        sum(sort(file_numbers)-(min(file_numbers):max(file_numbers))')~=0
    error(['There are MH draw files missing within chain ', int2str(blck)]) 
end


