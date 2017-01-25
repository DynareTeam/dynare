function acf = dyn_autocorr(y, ar)
% function acf = dyn_autocorr(y, ar)
% autocorrelation function of y
%
% INPUTS
%   y:       time series
%   ar:      # of lags
%
% OUTPUTS
%   acf:       autocorrelation for lags 1 to ar
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2015-16 Dynare Team
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


y=y(:);
acf = NaN(ar+1,1);
acf(1)=1;
m = mean(y);
sd = std(y,1);
for i=1:ar,
    acf(i+1) = (y(i+1:end)-m)'*(y(1:end-i)-m)./((size(y,1))*sd^2);
end
