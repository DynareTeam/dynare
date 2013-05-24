function plot_icforecast(Variables,periods,options_)
% Build plots for the conditional forecasts.
%
% INPUTS 
%  o Variables     [char]        m*x array holding the names of the endogenous variables to be plotted. 
%
% OUTPUTS
%  None.
% 
% SPECIAL REQUIREMENTS
%  This routine has to be called after imcforecast.m.

% Copyright (C) 2006-2012 Dynare Team
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

if exist('OCTAVE_VERSION') && octave_ver_less_than('3.4.0')
    % The set() command on the handle returned by area() crashes in Octave 3.2
    error('plot_conditional_forecast: you need Octave >= 3.4 (because of a bug in older versions)')
end

load conditional_forecasts;

if nargin==1 || isempty(periods) % Set default number of periods.
    eval(['periods = length(forecasts.cond.Mean.' Variables(1,:) ');']);
end

forecasts.graph.OutputDirectoryName = CheckPath('graphs',forecasts.graph.fname);

for i=1:size(Variables,1)
    eval(['ci1 = forecasts.cond.ci.' Variables(i,:) ';'])
    eval(['m1 = forecasts.cond.Mean.' Variables(i,:) ';'])
    eval(['ci2 = forecasts.uncond.ci.' Variables(i,:) ';'])
    eval(['m2 = forecasts.uncond.Mean.' Variables(i,:) ';'])
    build_figure(Variables(i,:),ci1(:,1:periods),ci2(:,1:periods),m1(1:periods),m2(1:periods),options_,forecasts.graph);
end

function build_figure(name,cci1,cci2,mm1,mm2,options_,graphoptions)
hh = dyn_figure(options_,'Name',['Conditional forecast (' graphoptions.title ,'): ' name '.']);
H = length(mm1);
h1 = area(1:H,cci1(2,1:H));
set(h1,'BaseValue',min([min(cci1(1,:)),min(cci2(1,:))]))
set(h1,'FaceColor',[.9 .9 .9])
hold on
h2 = area(1:H,cci1(1,1:H));
set(h2,'BaseValue',min([min(cci1(1,:)),min(cci2(1,:))]))
set(h2,'FaceColor',[1 1 1])
plot(1:H,mm1,'-k','linewidth',3)
plot(1:H,mm2,'--k','linewidth',3)
plot(1:H,cci2(1,:),'--k','linewidth',1)
plot(1:H,cci2(2,:),'--k','linewidth',1)
axis tight
hold off
dyn_saveas(hh,[graphoptions.OutputDirectoryName '/Conditional_forecast_',strrep(deblank(graphoptions.title),' ','_'),'_',deblank(name)],options_)