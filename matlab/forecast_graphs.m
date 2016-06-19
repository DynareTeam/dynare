function forecast_graphs(var_list,M_, oo_,options_)
% function forecast_graphs(var_list,M_, oo_,options_)
% Plots the classical forecasts created by dyn_forecast.m
% 
% Inputs:
%   o var_list              character array with variable names
%   o M_                    model structure
%   o oo_                   outputs structure
%   o options_              options structure

% Copyright (C) 2008-2016 Dynare Team
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

nc = 4;
nr = 3;
endo_names = M_.endo_names;
fname = M_.fname;
dname = M_.dname;

yf = oo_.forecast.Mean;
hpdinf = oo_.forecast.HPDinf;
hpdsup = oo_.forecast.HPDsup;
if isempty(var_list)
    var_list = endo_names(1:M_.orig_endo_nbr,:);
end
i_var = [];
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t an endogenous variable'])
    end
    i_var = [i_var; tmp];
end
nvar = length(i_var);

% create subdirectory <dname>/graphs if id doesn't exist
if ~exist(dname, 'dir')
    mkdir('.',dname);
end
if ~exist([dname '/graphs'],'dir')
    mkdir(dname,'graphs');
end

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX=write_LaTeX_header([M_.dname '/graphs/',fname,'_forcst.tex']);
end

m = 1;
n_fig = 1;
hh=dyn_figure(options_,'Name','Forecasts (I)');
for j= 1:nvar
    if m > nc*nr; 
        dyn_saveas(hh,[ dname '/graphs/forcst' int2str(n_fig)],options_);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s/graphs/forcst%d}\n',dname,n_fig);
            fprintf(fidTeX,'\\label{Fig:forcst:%d}\n',n_fig);
            fprintf(fidTeX,'\\caption{Mean forecasts and %2.0f\\%% confidence intervals}\n',options_.forecasts.conf_sig*100);
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
        end       
        n_fig =n_fig+1;
        eval(['hh=dyn_figure(options_,''Name'',''Forecasts (' int2str(n_fig) ')'');']);
        m = 1;
    end
    subplot(nr,nc,m);
    vn = deblank(endo_names(i_var(j),:));
    obs = 0;
% $$$         k = strmatch(vn,varobs,'exact'); 
% $$$   if ~isempty(k)
% $$$       yy = y.(vn)(end-9:end) + repmat(ys(i_var(j)),10,1)+trend(k,:)';
% $$$       plot(yy);
% $$$       hold on
% $$$       obs = 10;
% $$$   end
    plot([NaN(obs,1); yf.(vn)],'b-');
    hold on
    plot([NaN(obs,1); hpdinf.(vn)],'b-');
    hold on
    plot([NaN(obs,1); hpdsup.(vn)],'b-');
    title(vn,'Interpreter','none');
    xlim([1 obs+length(hpdsup.(vn))])
    hold off
    m = m + 1;
end

if m > 1
    dyn_saveas(hh,[dname '/graphs/forcst' int2str(n_fig)],options_);
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s/graphs/forcst%d}\n',dname,n_fig);
        fprintf(fidTeX,'\\label{Fig:forcst:%d}\n',n_fig);
        fprintf(fidTeX,'\\caption{Mean forecasts and %2.0f\\%% confidence intervals}\n',options_.forecasts.conf_sig*100);
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
end

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    write_LaTeX_footer(fidTeX);
end

if isfield(oo_.forecast,'HPDinf_ME')
    var_names=fieldnames(oo_.forecast.HPDinf_ME);
    
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX=write_LaTeX_header([M_.dname '/graphs/',fname,'_forcst_ME.tex']);
    end
    
    m = 1;
    n_fig = 1;
    hh=dyn_figure(options_,'Name','Forecasts including ME (I)');
    for j= 1:length(var_names)
        if m > nc*nr;
            dyn_saveas(hh,[ dname '/graphs/forcst_ME' int2str(n_fig)],options_);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s/graphs/forcst_ME%d}\n',options_.figures.textwidth*min((m-1)/nc,1),dname,n_fig);
                fprintf(fidTeX,'\\label{Fig:forcst_ME:%d}\n',n_fig);
                fprintf(fidTeX,'\\caption{Mean forecasts and %2.0f\\%% confidence intervals accounting for measurement error}\n',options_.conf_sig*100);
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
            end
            n_fig =n_fig+1;
            eval(['hh=dyn_figure(options_,''Name'',''Forecasts (' int2str(n_fig) ')'');']);
            m = 1;
        end
        subplot(nr,nc,m);
        vn = var_names{j};
        obs = 0;

        plot([NaN(obs,1); yf.(vn)],'b-');
        hold on
        plot([NaN(obs,1); oo_.forecast.HPDinf_ME.(vn)],'b-');
        hold on
        plot([NaN(obs,1); oo_.forecast.HPDsup_ME.(vn)],'b-');
        title(vn,'Interpreter','none');
        xlim([1 obs+length(oo_.forecast.HPDsup_ME.(vn))])
        hold off
        m = m + 1;
    end
    
    if m > 1
        dyn_saveas(hh,[dname '/graphs/forcst_ME' int2str(n_fig)],options_);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');          
            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s/graphs/forcst_ME%d}\n',options_.figures.textwidth*min((m-1)/nc,1),dname,n_fig);
            fprintf(fidTeX,'\\label{Fig:forcst_ME:%d}\n',n_fig);
            fprintf(fidTeX,'\\caption{Mean forecasts and %2.0f\\%% confidence intervals accounting for measurement error}\n',options_.conf_sig*100);
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
        end
    end
    
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        write_LaTeX_footer(fidTeX);
    end   
end


function fidTeX=write_LaTeX_header(fname)
%% write LaTeX-Header
fidTeX = fopen(fname,'w');
fprintf(fidTeX,'%% TeX eps-loader file generated by Dynare''s forecast_graphs.m.\n');
fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
fprintf(fidTeX,' \n');


function write_LaTeX_footer(fidTeX)
%% write LaTeX-Footer
fprintf(fidTeX,' \n');
fprintf(fidTeX,'%% End of TeX file.\n');
fclose(fidTeX);
