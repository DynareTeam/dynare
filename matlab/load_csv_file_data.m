function [freq, init, data, varlist] = load_csv_file_data(file, withtime, withnames, noemptycell)

%@info:
%! @deftypefn {Function File} {@var{freq}, @var{init}, @var{data}, @var{varlist} =} load_m_file_data (@var{file}, @var{withtime}, @var{withnames}, @var{noemptycell})
%! @anchor{load_csv_file_data}
%! @sp 1
%! Loads data in a csv file.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item file
%! string, name of the csv file
%! @item withtime
%! Scalar integer, non zero iff the first column is for the dates of the observations.
%! @item withnames
%! Scalar integer, non zero iff the first row is for the names of the variables.
%! @item noemptycell
%! Scalar integer, non zero iff the csv file does not have empty cells.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item freq
%! Scalar integer (1, 4, 12, 52).
%! @item init
%! dynDate object, initial date.
%! @item data
%! Matrix of doubles, data.
%! @item varlist
%! Cell of strings (names of the variables in the database).
%! @end table
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! The frequency and initial date can be specified with variables FREQ__ and INIT__ in the matlab/octave script. FREQ__ must be a scalar integer and INIT__ a string like '1938M11', '1945Q3', '1973W3' or '2009'. If these variables are not specified default values for freq and init are 1 and dynDate(1).
%! @end deftypefn
%@eod:

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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

% Set defaults.
if nargin<4
    noemptycell = 1;
    if nargin<3
        withnames = 1;
        if nargin <2
            withtime = 1;
            if nargin<1
                error('load_csv_file_data:: I need at least one input (name of the csv file)!')
            end
        end
    end
end

if ~withtime && ~withnames && noemptycell
    % Use matlab builtin routine!
    data = csvread(file);
end

if ~( isequal(withtime,0) || isequal(withtime,1) )
    error('load_csv_file_data:: Second input argument has to be equal to 1 or 0!')
end

if ~( isequal(withnames,0) || isequal(withnames,1) )
    error('load_csv_file_data:: Third input argument has to be equal to 1 or 0!')
end

% Output initialization
varlist = [];

% Check if file exists.
if check_file_extension(file,'csv')
    try
        fid = fopen(file,'r');
    catch
        error(['load_csv_file_data: I can''t find file ' file '!'])
    end
else
    error('load_csv_file_data: Wrong file extension!')
end

% bfile contains a vector of ascii codes.
bfile = fread(fid);

% Close (csv) file. 
fclose(fid);

% Set newline code (ok for *nix, check for mac and windows)
if isunix
    newline_code = 10;
elseif ispc
    newline_code = 13;
elseif ismac
    newline_code = 10;
else
    error('load_csv_file_data:: Not implemented for your OS!')
end

% Get the positions of the end-of-line code;
end_of_line_locations = find(bfile==newline_code);
if ispc && isempty(end_of_line_locations)
    newline_code=10;
    end_of_line_locations = find(bfile==newline_code);
end;
tmp = find(bfile==newline_code);

% Get the number of lines in the file.
ndx = length(tmp);                       

% Create a cell of indices for each line.
b = [1; end_of_line_locations+1];
c = [end_of_line_locations-1; length(bfile)+1];
b = b(1:end-1);
c = c(1:end-1);

linea = 1;

if withnames
    % Get the first line of the csv file (names of the variables).
    linee = char(transpose(bfile(b(linea):c(linea))));
    % Get the content of the first line and determine the number of variables and their names.
    [B,C] = get_cells_id(linee,',');
    if withtime
        B = B(2:end);
        C = C(2:end);
    end
    varlist = cell(length(B),1);
    number_of_variables = length(varlist);
    for i=1:number_of_variables
        varlist(i) = {linee(B(i):C(i))};
    end
    varlist = strtrim(varlist);
    linea = linea+1;
end

% Get following line (number 1 or 2 depending on withnames flag)
linee = char(transpose(bfile(b(linea):c(linea))));
comma_locations = transpose(strfind(linee,','));
B = 1;
C = comma_locations(1)-1;

if withtime
    tmp = linee(B:C);
    % Check the dates formatting
    if ~(isyearly(tmp) || isquaterly(tmp) || ismonthly(tmp) || isweekly(tmp))
        error('load_csv_file_data:: Formatting error. I can''t read the dates!')
    end
    if isyearly(tmp)==2
        % Remove the Y (gpm/iris date format) if necessary
        tmp = tmp(1:end-1);
    end
    init = dynDate(tmp);
    freq = init.freq;
    first = 2;
else
    init = dynDate(1);
    freq = 1;
    first = 1;
end

if ~withnames
    number_of_variables = length(tmp)-withtime;
end

% Initialization of matrix data.
data = zeros(ndx,number_of_variables);

% Populate data.
for linea = 1+withnames:ndx
    linee = char(transpose(bfile(b(linea):c(linea))));
    [B,C] = get_cells_id(linee,',');
    for i=first:length(B)
        data(linea,i-withtime) = str2double(linee(B(i):C(i)));
    end
end

% Remove first line if withnames
data = data(1+withnames:ndx,:);