function trace_plot(options_,M_,estim_params_,type,blck,name1,name2)
% This function builds trace plot for the Metropolis-Hastings draws.
%
% INPUTS
%
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   estim_params_   [structure]    Dynare structure (related to estimation).
%   type            [string]       'DeepParameter', 'MeasurementError' (for measurement equation error),
%                                  'StructuralShock' (for structural shock)
%                                  or 'PosteriorDensity (for posterior density)'
%   blck            [integer]      Number of the mh chain.
%   name1           [string]       Object name.
%   name2           [string]       Object name.
%
% OUTPUTS
%   None
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2016 Dynare Team
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

% Cet the column index:
if strcmpi(type,'PosteriorDensity')
    column=0;
    name1='';
else
    if nargin<7
        column = name2index(options_, M_, estim_params_, type, name1);
    else
        column = name2index(options_, M_, estim_params_, type, name1, name2);
    end
end

if isempty(column)
    return
end

% Get informations about the posterior draws:
MetropolisFolder = CheckPath('metropolis',M_.dname);
load_last_mh_history_file(MetropolisFolder, M_.fname);

FirstMhFile = 1;
FirstLine = 1;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
clear record;

% Get all the posterior draws:
PosteriorDraws = GetAllPosteriorDraws(column, FirstMhFile, FirstLine, TotalNumberOfMhFiles, TotalNumberOfMhDraws, blck);


% Plot the posterior draws:

if strcmpi(type,'DeepParameter')
    TYPE = 'parameter ';
elseif strcmpi(type,'StructuralShock')
    if nargin<7
        TYPE = 'the standard deviation of structural shock ';
    else
        TYPE = 'the correlation between structural shocks ';
    end
elseif strcmpi(type,'MeasurementError')
    if nargin<7
        TYPE = 'the standard deviation of measurement error ';
    else
        TYPE = 'the correlation between measurement errors ';
    end
elseif strcmpi(type,'PosteriorDensity')
    TYPE='the posterior density';
end

if nargin<7
    FigureName = ['Trace plot for ' TYPE name1];
else
    FigureName = ['Trace plot for ' TYPE name1 ' and ' name2];
end

if options_.mh_nblck>1
    FigureName = [ FigureName , ' (block number ' int2str(blck)  ').'];
end

hh=dyn_figure(options_,'Name',FigureName);
plot(1:TotalNumberOfMhDraws,PosteriorDraws,'Color',[.7 .7 .7]);


% Compute the moving average of the posterior draws:

N = options_.trace_plot_ma;
MovingAverage = NaN(TotalNumberOfMhDraws,1);
first = N+1;
last = TotalNumberOfMhDraws-N;

for t=first:last
    MovingAverage(t) = mean(PosteriorDraws(t-N:t+N));
end

hold on
plot(1:TotalNumberOfMhDraws,MovingAverage,'-k','linewidth',2)
hold off
axis tight
legend({'MCMC draw';[num2str(N) ' period moving average']},'Location','NorthWest')
% create subdirectory <fname>/graphs if it doesn't exist
if ~exist(M_.fname, 'dir')
    mkdir('.',M_.fname);
end
if ~exist([M_.fname filesep 'graphs'],'dir')
    mkdir(M_.fname,'graphs');
end

%get name for plot
if strcmpi(type,'PosteriorDensity')
    plot_name='Posterior';
else
    plot_name=get_the_name(column,0,M_,estim_params_,options_);
end
plot_name=[plot_name,'_blck_',num2str(blck)];

dyn_saveas(hh,[M_.fname, filesep, 'graphs', filesep, 'TracePlot_' plot_name],options_)

if options_.TeX
    fid=fopen([M_.fname,'/graphs/',M_.fname,'_TracePlot_' plot_name,'.tex'],'w+');
    
    if strcmpi(type,'DeepParameter')
        tex_names=M_.param_names_tex;
        base_names=M_.param_names;
    elseif strcmpi(type,'StructuralShock')
        tex_names=M_.exo_names_tex;
        base_names=M_.exo_names;
    elseif strcmpi(type,'MeasurementError')
        tex_names=M_.endo_names_tex;
        base_names=M_.endo_names;
    end
    
    if strcmpi(type,'PosteriorDensity')
        FigureName = ['Trace plot for ' TYPE name1];
    else
        if nargin<7
            FigureName = ['Trace plot for ' TYPE '$' deblank(tex_names(strmatch(name1,base_names,'exact'),:)) '$'];
        else
            FigureName = ['Trace plot for ' TYPE '$' deblank(tex_names(strmatch(name1,base_names,'exact'),:)) '$ and $' deblank(tex_names(strmatch(name2,base_names,'exact'),:)) '$'];
        end
    end
    if options_.mh_nblck>1
        FigureName = [ FigureName , ' (block number ' int2str(blck)  ').'];
    end
    
    fprintf(fid,'%-s\n','\begin{figure}[H]');
    fprintf(fid,'%-s\n','\centering');
    fprintf(fid,'%-s\n',['  \includegraphics[width=0.8\textwidth]{',[M_.fname, '/graphs/TracePlot_' plot_name],'}\\']);
    fprintf(fid,'%-s\n',['    \caption{',FigureName,'}']);
    fprintf(fid,'%-s\n','\end{figure}');
    fclose(fid);
end

