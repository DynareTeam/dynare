function list_of_exported_variables_ = dat_fil_(dat_fil_to_load_);
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.

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

try
  eval(dat_fil_to_load_);
catch
  load(dat_fil_to_load_);
end
clear dat_fil_to_load_;

list_of_local_variables_=who;

for j=1:length(list_of_local_variables_),
  eval(['list_of_exported_variables_.',list_of_local_variables_{j},'=',list_of_local_variables_{j},';']);
end