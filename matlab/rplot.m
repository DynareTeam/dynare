function rplot(s1)
% function rplot(s1)
%
% Plots the simulated trajectory of one or several variables.
% The entire simulation period is plotted, unless instructed otherwise
% with "dsample".
%
% INPUTS
%    s1:           character matrix of variable names
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2016 Dynare Team
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

global M_ oo_ options_

rplottype = options_.rplottype;

if isempty(oo_.endo_simul)
    error('rplot: oo_.endo_simul is empty.')
end

% create subdirectory <fname>/graphs if it doesn't exist
if ~exist(M_.fname, 'dir')
    mkdir('.',M_.fname);
end
if ~exist([M_.fname filesep 'graphs'],'dir')
    mkdir(M_.fname,'graphs');
end

col = ['y','c','r','g','b','w','m'] ;
ix = [1 - M_.maximum_lag:size(oo_.endo_simul,2)-M_.maximum_lag]' ;

y = [];
for k=1:size(s1,1)
    if isempty(strmatch(deblank(s1(k,:)),M_.endo_names,'exact')) 
        if isempty(strmatch(deblank(s1(k,:)),M_.exo_names,'exact')) 
            error (['rplot: One of the variables specified does not exist']) ;
        else
            y = [y; oo_.exo_simul(:,strmatch(deblank(s1(k,:)),M_.exo_names,'exact'))'] ;        
        end
    else
        y = [y; oo_.endo_simul(strmatch(deblank(s1(k,:)),M_.endo_names,'exact'),:)] ;
    end
end

if options_.smpl == 0
    i = [max(1, M_.maximum_lag):size(oo_.endo_simul,2)]' ;
else
    i = [options_.smpl(1)+M_.maximum_lag:options_.smpl(2)+M_.maximum_lag]' ;
end

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([M_.fname, filesep, 'graphs', filesep, M_.fname '_simulated_trajectories_', num2str(rplottype), '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by rplot.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
end


t = ['Plot of '] ;
if rplottype == 0
    for j = 1:size(y,1)
        t = [t s1(j,:) ' '] ;
    end
    hh=dyn_figure(options_,'Name',['Simulated Trajectory']);
    plot(ix(i),y(:,i)) ;
    title (t,'Interpreter','none') ;
    xlabel('Periods') ;
    if size(s1,1) > 1
        if isoctave
            legend(s1);
        else
            h = legend(s1);
            set(h, 'Interpreter', 'none');
        end
    end
    dyn_saveas(hh,[M_.fname, filesep, 'graphs', filesep, 'SimulatedTrajectory_' deblank(s1(1,:))],options_)
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        create_TeX_loader(fidTeX,options_,[M_.fname, '/graphs/', 'SimulatedTrajectory_' deblank(s1(1,:))],'Simulated trajectories','SimulatedTrajectory_',deblank(s1(1,:)),1)
    end
elseif rplottype == 1
    for j = 1:size(y,1)
        hh=dyn_figure(options_,'Name',['Simulated Trajectory']);
        plot(ix(i),y(j,i)) ;
        title(['Plot of ' s1(j,:)],'Interpreter','none') ;
        xlabel('Periods') ;
        dyn_saveas(hh,[M_.fname, filesep, 'graphs', filesep, 'SimulatedTrajectory_' deblank(s1(j,:))],options_)
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            create_TeX_loader(fidTeX,options_,[M_.fname, '/graphs/', 'SimulatedTrajectory_' deblank(s1(j,:))],'Simulated trajectories','SimulatedTrajectory_',deblank(s1(j,:)),1);
        end
    end
elseif rplottype == 2
    hh=dyn_figure(options_,'Name',['Simulated Trajectory']);
    nl = max(1,fix(size(y,1)/4)) ;
    nc = ceil(size(y,1)/nl) ;
    for j = 1:size(y,1)
        subplot(nl,nc,j) ;
        plot(ix(i),y(j,i)) ;
        hold on ;
        if ~isempty(strmatch(deblank(s1(j,:)),M_.endo_names,'exact'))
            plot(ix(i),oo_.steady_state(strmatch(deblank(s1(j,:)),M_.endo_names,'exact'))*ones(1,size(i,1)),'r:') ;
        else
            plot(ix(i),oo_.exo_steady_state(strmatch(deblank(s1(j,:)),M_.exo_names,'exact'))*ones(1,size(i,1)),'r:') ;
        end
        xlabel('Periods') ;
        ylabel([s1(j,:)],'Interpreter','none') ;
        title(['Plot of ' s1(j,:)],'Interpreter','none') ;
        axis tight;
    end
    dyn_saveas(hh,[M_.fname, filesep, 'graphs', filesep, 'SimulatedTrajectory_' deblank(s1(1,:))],options_)
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        create_TeX_loader(fidTeX,options_,[M_.fname, '/graphs/', 'SimulatedTrajectory_' deblank(s1(1,:))],'Simulated trajectories','SimulatedTrajectory_',deblank(s1(1,:)),min(j/nc,1));
    end
end

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end

function []=create_TeX_loader(fidTeX,options,figpath,caption,label_name,label_type,scale_factor)
    if nargin<6
        scale_factor=1;
    end
    fprintf(fidTeX,' \n'); 
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',0.8*scale_factor,strrep(figpath,'\','/'));
    fprintf(fidTeX,'\\caption{%s.}',caption);
    fprintf(fidTeX,'\\label{Fig:%s:%s}\n',label_name,label_type);
    fprintf(fidTeX,'\\end{figure}\n\n');

% 02/28/01 MJ replaced bseastr by MATLAB's strmatch
% 06/19/01 MJ added 'exact' to strmatch calls
% 06/25/03 MJ correction when options_.smpl ~= 0
% 03/18/13 JP bugfix for rplottype>0; added figure names