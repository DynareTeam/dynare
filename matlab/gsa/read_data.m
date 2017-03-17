function [gend, data] = read_data()
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2012-2015 European Commission
% Copyright (C) 2012-2015 Dynare Team
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

global options_

rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);

options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;

rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
if options_.loglinear == 1 & ~options_.logdata
  rawdata = log(rawdata);
end
if options_.prefilter == 1
  data = transpose(rawdata-ones(gend,1)* mean(rawdata,1));
else
  data = transpose(rawdata);
end

if ~isreal(rawdata)
  error(['There are complex values in the data. Probably  a wrong' ...
	 ' transformation'])
end
