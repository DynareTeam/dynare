function [freq,init,data,varlist,tex] = load_mat_file_data(file)  % --*-- Unitary tests --*--

% Loads data in a matlab/octave mat-file.
%
% INPUTS 
%  o file         string, name of the matlab/octave mat file (with path)
%
% OUTPUTS 
%  o freq        integer scalar equal to 1, 4, 12 or 52 (for annual, quaterly, monthly or weekly frequencies).
%  o init        dates object, initial date in the dataset.
%  o data        matrix of doubles, the data.
%  o varlist     cell of strings, names of the variables.
%
% REMARKS 
% The frequency and initial date can be specified with variables FREQ__ and INIT__ in the matlab/octave binary file. FREQ__ must 
% be a scalar integer and INIT__ a string like '1938M11', '1945Q3', '1973W3' or '2009A'. If these variables are not specified 
% default values for freq and init are 1 and dates(1,1).

% Copyright (C) 2012-2014 Dynare Team
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

datafile = load(file);

if isfield(datafile,'INIT__')
    if isdate(datafile.INIT__)
        init = dates(datafile.INIT__);
        datafile = rmfield(datafile, 'INIT__');
    else
        error('load_mat_file_data: INIT__ cannot be interpreted as a date.')
    end
else
    init = dates(1,1);
end

if isfield(datafile,'FREQ__')
    freq = datafile.FREQ__;
    datafile = rmfield(datafile, 'FREQ__');
else
    freq = init.freq;
end

if ~isequal(freq,init.freq)
    error('load_mat_file_data: INIT__ and FREQ__ are not consistent!')
end

if isfield(datafile,'NAMES__')
    varlist = datafile.NAMES__;
    datafile = rmfield(datafile, 'NAMES__');
else
    varlist = [];
end

if isfield(datafile,'TEX__')
    tex = datafile.TEX__;
    datafile = rmfield(datafile, 'TEX__');
else
    tex = [];
end

data = [];
if isempty(varlist)
    varlist = fieldnames(datafile);
end

for i=1:length(varlist)
    try
        tmp = getfield(datafile,varlist{i});
        data = [data,  tmp(:)];
    catch
        error(['load_mat_file:: All the vectors (variables) in ' inputname(1) ' must have the same number of rows (observations)!'])
    end
end

%@test:1
%$ % Create a data mat-file
%$ FREQ__ = 12;
%$ INIT__ = '1938M11';
%$ NAMES__ = {'hagop'; 'bedros'};
%$ TEX__ = NAMES__;
%$ hagop  = [1; 2; 3; 4; 5];
%$ bedros = [2; 3; 4; 5; 6];
%$ save('datafile_for_test');
%$
%$ % Try to read the data mat-file
%$ t = zeros(8,1);
%$ try
%$     [freq,init,data,varlist,tex] = load_mat_file_data('datafile_for_test');
%$     t(1) = 1;
%$ catch exception
%$     t = t(1);
%$     T = all(t);
%$     LOG = getReport(exception,'extended');
%$     return
%$ end
%$
%$ % Check the results.
%$ t(2) = dassert(freq,12);
%$ t(3) = dassert(isa(init,'dates'),true);
%$ t(4) = dassert(init.freq,12);
%$ t(5) = dassert(init.time,[1938 11]);
%$ t(6) = dassert(varlist,{'hagop';'bedros'});
%$ t(7) = dassert(varlist,{'hagop';'bedros'});
%$ t(8) = dassert(data(:,1),[1;2;3;4;5]);
%$ t(9) = dassert(data(:,2),[2;3;4;5;6]);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Create a data mat-file
%$ FREQ__ = 12;
%$ INIT__ = '1938M11';
%$ NAMES__ = {'hagop'; 'bedros'};
%$ TEX__ = NAMES__;
%$ hagop  = [1, 2, 3, 4, 5];
%$ bedros = [2, 3, 4, 5, 6];
%$ save('datafile_for_test');
%$
%$ % Try to read the data mat-file
%$ t = zeros(8,1);
%$ try
%$     [freq,init,data,varlist,tex] = load_mat_file_data('datafile_for_test');
%$     t(1) = 1;
%$ catch exception
%$     t = t(1);
%$     T = all(t);
%$     LOG = getReport(exception,'extended');
%$     return
%$ end
%$
%$ % Check the results.
%$ t(2) = dassert(freq,12);
%$ t(3) = dassert(isa(init,'dates'),true);
%$ t(4) = dassert(init.freq,12);
%$ t(5) = dassert(init.time,[1938 11]);
%$ t(6) = dassert(varlist,{'hagop';'bedros'});
%$ t(7) = dassert(varlist,{'hagop';'bedros'});
%$ t(8) = dassert(data(:,1),[1;2;3;4;5]);
%$ t(9) = dassert(data(:,2),[2;3;4;5;6]);
%$ T = all(t);
%@eof:2
