function sout = myboxplot (data,notched,symbol,vertical,maxwhisker)
%  sout = myboxplot (data,notched,symbol,vertical,maxwhisker)

% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2012 European Commission
% Copyright (C) 2012 Dynare Team
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

% % % % endif
if nargin < 5 | isempty(maxwhisker), maxwhisker = 1.5; end
if nargin < 4 | isempty(vertical), vertical = 1; end
if nargin < 3 | isempty(symbol), symbol = ['+','o']; end
if nargin < 2 | isempty(notched), notched = 0; end

if length(symbol)==1, symbol(2)=symbol(1); end

if notched==1, notched=0.25; end
a=1-notched;

% ## figure out how many data sets we have
if iscell(data), 
  nc = length(data);
else
%   if isvector(data), data = data(:); end
  nc = size(data,2);
end

% ## compute statistics
% ## s will contain
% ##    1,5    min and max
% ##    2,3,4  1st, 2nd and 3rd quartile
% ##    6,7    lower and upper confidence intervals for median
s = zeros(7,nc);
box = zeros(1,nc);
whisker_x = ones(2,1)*[1:nc,1:nc];
whisker_y = zeros(2,2*nc);
outliers_x = [];
outliers_y = [];
outliers2_x = [];
outliers2_y = [];

for i=1:nc
    %   ## Get the next data set from the array or cell array
    if iscell(data)
        col = data{i}(:);
    else
        col = data(:,i);
    end
%   ## Skip missing data
% % % % % % %   col(isnan(col) | isna (col)) = [];
  col(isnan(col)) = [];

  %   ## Remember the data length
  nd = length(col);
  box(i) = nd;
  if (nd > 1)
%     ## min,max and quartiles
%     s(1:5,i) = statistics(col)(1:5);
s(1,i)=min(col);
s(5,i)=max(col);
s(2,i)=myprctilecol(col,25);
s(3,i)=myprctilecol(col,50);
s(4,i)=myprctilecol(col,75);








%     ## confidence interval for the median
    est = 1.57*(s(4,i)-s(2,i))/sqrt(nd);
    s(6,i) = max([s(3,i)-est, s(2,i)]);
    s(7,i) = min([s(3,i)+est, s(4,i)]);
%     ## whiskers out to the last point within the desired inter-quartile range
    IQR = maxwhisker*(s(4,i)-s(2,i));
    whisker_y(:,i) = [min(col(col >= s(2,i)-IQR)); s(2,i)];
    whisker_y(:,nc+i) = [max(col(col <= s(4,i)+IQR)); s(4,i)];
%     ## outliers beyond 1 and 2 inter-quartile ranges
    outliers = col((col < s(2,i)-IQR & col >= s(2,i)-2*IQR) | (col > s(4,i)+IQR & col <= s(4,i)+2*IQR));
    outliers2 = col(col < s(2,i)-2*IQR | col > s(4,i)+2*IQR);
    outliers_x = [outliers_x; i*ones(size(outliers))];
    outliers_y = [outliers_y; outliers];
    outliers2_x = [outliers2_x; i*ones(size(outliers2))];
    outliers2_y = [outliers2_y; outliers2];
  elseif (nd == 1)
%     ## all statistics collapse to the value of the point
    s(:,i) = col;
%     ## single point data sets are plotted as outliers.
    outliers_x = [outliers_x; i];
    outliers_y = [outliers_y; col];
  else
%     ## no statistics if no points
    s(:,i) = NaN;
  end
end
% % % % if isempty(outliers2_y)
% % % %     outliers2_y=
% ## Note which boxes don't have enough stats
chop = find(box <= 1);
    
% ## Draw a box around the quartiles, with width proportional to the number of
% ## items in the box. Draw notches if desired.
box = box*0.23/max(box);
quartile_x = ones(11,1)*[1:nc] + [-a;-1;-1;1;1;a;1;1;-1;-1;-a]*box;
quartile_y = s([3,7,4,4,7,3,6,2,2,6,3],:);

% ## Draw a line through the median
median_x = ones(2,1)*[1:nc] + [-a;+a]*box;
% median_x=median(col);
median_y = s([3,3],:);

% ## Chop all boxes which don't have enough stats
quartile_x(:,chop) = [];
quartile_y(:,chop) = [];
whisker_x(:,[chop,chop+nc]) = [];
whisker_y(:,[chop,chop+nc]) = [];
median_x(:,chop) = [];
median_y(:,chop) = [];
% % % % 
% ## Add caps to the remaining whiskers
cap_x = whisker_x;
cap_x(1,:) =cap_x(1,:)- 0.05;
cap_x(2,:) =cap_x(2,:)+ 0.05;
cap_y = whisker_y([1,1],:);

% #quartile_x,quartile_y
% #whisker_x,whisker_y
% #median_x,median_y
% #cap_x,cap_y
% 
% ## Do the plot

mm=min(min(data));
MM=max(max(data));

if vertical
    plot (quartile_x, quartile_y, 'b',  ...
        whisker_x, whisker_y, 'b--',   ...
        cap_x, cap_y, 'k',   ...
        median_x, median_y, 'r',  ...
        outliers_x, outliers_y, [symbol(1),'r'],   ...
        outliers2_x, outliers2_y, [symbol(2),'r']);
        set(gca,'XTick',1:nc);
        set(gca, 'XLim', [0.5, nc+0.5]);
        set(gca, 'YLim', [mm-(MM-mm)*0.05-eps, MM+(MM-mm)*0.05+eps]);

else
% % % % %     plot (quartile_y, quartile_x, "b;;",
% % % % %     whisker_y, whisker_x, "b;;",
% % % % %     cap_y, cap_x, "b;;",
% % % % %     median_y, median_x, "r;;",
% % % % %     outliers_y, outliers_x, [symbol(1),"r;;"],
% % % % %     outliers2_y, outliers2_x, [symbol(2),"r;;"]);
end

if nargout,
  sout=s;
end
% % % endfunction