function dyn_graph=dynare_graph(dyn_graph,y,tit,x)
% function dynare_graph(y,tit,x) 
% graphs
%
% INPUT
%   figure_name: name of the figures
%   colors: line colors
%
% OUTPUT
%   dyn_graph: structure with figure information
%
% SPECIAL REQUIREMENT
%   none

% Copyright (C) 2006-2013 Dynare Team
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

if nargin < 4
    x = (1:size(y,1))';
end
nplot = dyn_graph.plot_nbr + 1; 
if nplot > dyn_graph.max_nplot
    figure('Name',dyn_graph.figure_name);
    nplot = 1;
end
dyn_graph.plot_nbr = nplot;
subplot(dyn_graph.nr,dyn_graph.nc,nplot);

line_types = dyn_graph.line_types;
line_type = line_types{1};
for i=1:size(y,2);
    if length(line_types) > 1
        line_type = line_types{i};
    end
    
    plot(x,y(:,i),line_type);
    hold on
end
title(tit);
hold off