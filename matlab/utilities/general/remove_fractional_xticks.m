function remove_fractional_xticks
% function remove_fractional_xticks
% removes non-integer xtick-labels

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

xticks=get(gca,'xtick');
fractional_periods=find(rem(xticks,1)~=0);
if ~isempty(fractional_periods)
    xticks(fractional_periods)=[];
    xticklabels=get(gca,'xticklabel');
    xticklabels(fractional_periods)=[];
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',xticklabels);
end