function dyn_data_01=read_variables(file_name_01,var_names_01,dyn_data_01,xls_sheet,xls_range)

% function dyn_data_01=read_variables(file_name_01,var_names_01,dyn_data_01,xls_sheet,xls_range)
% Read data
%
% INPUTS
%    file_name_01:    file name
%    var_names_01:    variables name
%    dyn_data_01:     
%    xls_sheet:       Excel sheet name
%    xls_range:       Excel range specification
%
% OUTPUTS
%    dyn_data_01:
%
% SPECIAL REQUIREMENTS
% all local variables have complicated names in order to avoid name
% conflicts with possible user variable names

% Copyright (C) 2005-2013 Dynare Team
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


old_pwd = pwd;
[directory,basename,extension] = fileparts(file_name_01);
if ~isempty(directory)
    cd(directory)
end

dyn_size_01 = size(dyn_data_01,1);
var_size_01 = size(var_names_01,1);

% Auto-detect extension if not provided
if isempty(extension)
    if exist([basename '.m'],'file')
        extension = '.m';
    elseif exist([basename '.mat'],'file')
        extension = '.mat';
    elseif exist([basename '.xls'],'file')
        extension = '.xls';
    elseif exist([basename '.xlsx'],'file')
        extension = '.xlsx';
    elseif exist([basename '.csv'],'file')
        extension = '.csv';
    else
        error(['Can''t find datafile: ' basename '.{m,mat,xls,xlsx,csv}']);
    end
end

fullname = [basename extension];

if ~exist(fullname)
    error(['Can''t find datafile: ' fullname ]);
end 

switch (extension)
    case '.m'
        eval(basename);
        for dyn_i_01=1:var_size_01
            dyn_tmp_01 = eval(var_names_01(dyn_i_01,:));
            if length(dyn_tmp_01) > dyn_size_01 && dyn_size_01 > 0
                cd(old_pwd)
                error('data size is too large')
            end
            dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
        end
    case '.mat'
        s = load(basename);
        for dyn_i_01=1:var_size_01
            dyn_tmp_01 = s.(deblank(var_names_01(dyn_i_01,:)));
            if length(dyn_tmp_01) > dyn_size_01 && dyn_size_01 > 0
                cd(old_pwd)
                error('data size is too large')
            end
            dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
        end
    case { '.xls', '.xlsx' }
        [num,txt,raw] = xlsread(fullname,xls_sheet,xls_range); % Octave needs the extension explicitly
        for i=1:size(raw,2)
            if isnan(raw{1,i})
                raw{1,i} = ' ';
            end
        end
        for dyn_i_01=1:var_size_01
            iv = strmatch(var_names_01(dyn_i_01,:),raw(1,:),'exact');
            dyn_tmp_01 = [raw{2:end,iv}]';
            if length(dyn_tmp_01) > dyn_size_01 && dyn_size_01 > 0
                cd(old_pwd)
                error('data size is too large')
            end
            dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
        end
    case '.csv'
        [freq,init,data,varlist] = load_csv_file_data(fullname);
        %var_names_01 = deblank(var_names_01);
        for dyn_i_01=1:var_size_01
            iv = strmatch(strtrim(var_names_01(dyn_i_01,:)),varlist,'exact');
            if ~isempty(iv)
                dyn_tmp_01 = [data(2:end,iv)]';
                if length(dyn_tmp_01) > dyn_size_01 && dyn_size_01 > 0
                    cd(old_pwd)
                    error('data size is too large')
                end
                dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
            else
                cd(old_pwd)
                error([strtrim(var_names_01(dyn_i_01,:)) ' not found in ' fullname])
            end
        end
    otherwise
        cd(old_pwd)
        error(['Unsupported extension for datafile: ' extension])
end

cd(old_pwd)
disp(sprintf('Loading %d observations from %s\n',...
             size(dyn_data_01,1),fullname))
