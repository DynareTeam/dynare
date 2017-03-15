function expand_group(use_shock_groups,var_list_, ic)
% function expand_group(use_shock_groups,var_list_, ic)
% Expands shocks contributions out of a group of shocks
%
% INPUTS
%    use_shock_groups: [char]  name of the group
%    var_list_:        [char]  list of variables
%    ic:               [int]   group # to expand
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2016-2017 Dynare Team
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

M = evalin('base','M_');
oo = evalin('base','oo_');
options = evalin('base','options_');
mydata=get(findobj(gcf,'tag',['group' int2str(ic)]),'userdata');
if isfield(mydata,'shock_decomp')
    options.shock_decomp=mydata.shock_decomp;
end
% define expanded group
label = mydata.shock_group.label;
shocks = mydata.shock_group.shocks;
options.shock_decomp.fig_names = [mydata.fig_names '. Expand'];
options.use_shock_groups = strrep(label,' ','_'); %[use_shock_groups_old int2str(ic)];
for j=1:length(shocks)
    M.shock_groups.(options.use_shock_groups).(['group' int2str(j)]).label=shocks{j};
    M.shock_groups.(options.use_shock_groups).(['group' int2str(j)]).shocks=shocks(j);
end 

options.shock_decomp.interactive=0;
options.shock_decomp.expand=1;
plot_shock_decomposition(M,oo,options,var_list_);


