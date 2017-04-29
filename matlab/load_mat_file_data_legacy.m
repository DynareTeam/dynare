function data  = load_mat_file_data_legacy(datafile, varobs)

% Copyright (C) 2017 Dynare Team
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

data_file=load(datafile);

names=fieldnames(data_file);

if ~all(ismember(varobs',names))
    missing_variables=varobs(~ismember(varobs',names))';
    disp_string=[missing_variables{1,:}];
    for ii=2:size(missing_variables,1)
        disp_string=[disp_string,', ',missing_variables{ii,:}];
    end
    error('makedataset: The variable(s) %s listed in varobs are not contained in the dataset %s',disp_string);    
else  
    data_mat=[];
    for var_iter=1:length(varobs)
        try
           data_mat=[data_mat data_file.(varobs{1,var_iter})];
        catch
            error('makedataset: The variable %s does not have dimensions conformable with the previous one',varobs{1,var_iter});                
        end
    end
end
data = dseries(data_mat,[],varobs);