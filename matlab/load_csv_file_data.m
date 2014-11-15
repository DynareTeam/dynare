function [freq, init, data, varlist] = load_csv_file_data(file)

% Loads data in a csv file.
%
% INPUTS 
%  o file        string, name of the csv file (with path).
%
% OUTPUTS 
%  o freq        integer scalar equal to 1, 4, 12 or 52 (for annual, quaterly, monthly or weekly frequencies).
%  o init        dates object, initial date in the dataset.
%  o data        matrix of doubles, the data.
%  o varlist     cell of strings, names of the variables.
%
% REMARKS
%  The varlist output will be set only if the first line contains variable
%  names. Similarly, if the first column does not contain dates, then
%  freq will be 1 and init will be year 1.

% Copyright (C) 2012-2013 Dynare Team
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

% Output initialization
freq = 1;                  % Default frequency is annual.
init = dates(1,1);         % Default initial date is year one.
varlist = [];

assert(exist(file, 'file') == 2, ['load_csv_file_data: I can''t find file ' file '!']);

if isoctave
    if ~user_has_octave_forge_package('io')
        error('The io package is required to read CSV files from Octave')
    end
    A = csv2cell(file);
    [data, T, L] = parsecell(A);
    withvars = L.numlimits(2,1) > L.txtlimits(2,1);
    withtime = L.numlimits(1,1) > L.txtlimits(1,1);
else
    A = importdata(file, ',');
    if ~isstruct(A)
        data = A;
        T = {};
        withvars = 0;
        withtime = 0;
    else
        data = A.data;
        T = A.textdata;
        % importdata() allows text only at the top and the left, so the following
        %  tests are sufficient.
        withvars = size(T, 2) >= size(data, 2);
        withtime = size(T, 1) >= size(data, 1);
    end
end

if withvars
    varlist = T(1, 2:end);
    T = T(2:end, :);
end
if withtime
    init = dates(T{1, 1});
    freq = init.freq;
end

varlist = transpose(varlist);
