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

% Copyright (C) 2001-2009 Dynare Team
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

col = ['y','c','r','g','b','w','m'] ;
ix = [1 - M_.maximum_lag:size(oo_.endo_simul,2)-M_.maximum_lag]' ;

y = [];
for k=1:size(s1,1)
    if isempty(strmatch(s1(k,:),M_.endo_names,'exact'))
        error (['One of the variable specified does not exist']) ;
    end

    y = [y; oo_.endo_simul(strmatch(s1(k,:),M_.endo_names,'exact'),:)] ;
end

if options_.smpl == 0
    i = [max(1, M_.maximum_lag):size(oo_.endo_simul,2)]' ;
else
    i = [options_.smpl(1)+M_.maximum_lag:options_.smpl(2)+M_.maximum_lag]' ;
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
        if exist('OCTAVE_VERSION')
            legend(s1, 0);
        else
            h = legend(s1,0);
            set(h, 'Interpreter', 'none');
        end
    end
elseif rplottype == 1
    for j = 1:size(y,1)
        hh=dyn_figure(options_,'Name',['Simulated Trajectory']);
        plot(ix(i),y(j,i)) ;
        title(['Plot of ' s1(j,:)],'Interpreter','none') ;
        xlabel('Periods') ;
    end
elseif rplottype == 2
    hh=dyn_figure(options_,'Name',['Simulated Trajectory']);
    nl = max(1,fix(size(y,1)/4)) ;
    nc = ceil(size(y,1)/nl) ;
    for j = 1:size(y,1)
        subplot(nl,nc,j) ;
        plot(ix(i),y(j,i)) ;
        hold on ;
        plot(ix(i),oo_.steady_state(strmatch(s1(j,:),M_.endo_names,'exact'))*ones(1,size(i,1)),'r:') ;
        xlabel('Periods') ;
        ylabel([s1(j,:)],'Interpreter','none') ;
        title(['Plot of ' s1(j,:)],'Interpreter','none') ;
        axis tight;
    end
end
% create subdirectory <fname>/graphs if it doesn't exist
if ~exist(M_.fname, 'dir')
    mkdir('.',M_.fname);
end
if ~exist([M_.fname filesep 'graphs'])
    mkdir(M_.fname,'graphs');
end

dyn_saveas(hh,[M_.fname, filesep, 'graphs', filesep, 'SimulatedTrajectory_' deblank(s1(1,:))],options_)


% 02/28/01 MJ replaced bseastr by MATLAB's strmatch
% 06/19/01 MJ added 'exact' to strmatch calls
% 06/25/03 MJ correction when options_.smpl ~= 0
% 03/18/13 JP bugfix for rplottype>0; added figure names









