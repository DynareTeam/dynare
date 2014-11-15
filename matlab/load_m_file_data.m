function [freq,init,data,varlist,tex] = load_m_file_data(file)

% Loads data in a matlab/octave script.
%
% INPUTS 
%  o file         string, name of the matlab/octave script (with path)
%
% OUTPUTS 
%  o freq        integer scalar equal to 1, 4, 12 or 52 (for annual, quaterly, monthly or weekly frequencies).
%  o init        dates object, initial date in the dataset.
%  o data        matrix of doubles, the data.
%  o varlist     cell of strings, names of the variables.
%
% REMARKS 
% The frequency and initial date can be specified with variables FREQ__ and INIT__ in the matlab/octave script. FREQ__ must 
% be a scalar integer and INIT__ a string like '1938M11', '1945Q3', '1973W3' or '2009A'. If these variables are not specified 
% default values for freq and init are 1 and dates(1,1).

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

if isoctave
    run(file);
else
    basename = file(1:end-2);
    run(basename);
end

if exist('INIT__','var')
    if isdate(INIT__)
        init = dates(INIT__);
        clear('INIT__')
    else
        error('load_m_file_data: INIT__ cannot be interpreted as a date.')
    end
else
    init = dates(1,1); % Default initial date is year one.
end

if exist('FREQ__','var')
    freq = FREQ__;
    clear('FREQ__');
else
    freq = init.freq;
end

if ~isequal(freq,init.freq)
    error('load_m_file_data: INIT__ and FREQ__ are not consistent!')
end

if exist('NAMES__','var')
    varlist0 = NAMES__;
    clear('NAMES__');
else
    varlist0 = [];
    list_of_variables = [];
end

if exist('TEX__','var')
    tex = TEX__;
    clear('TEX__');
else
    tex = [];
end


if isempty(varlist0)
    list_of_variables = whos();
end

data = [];
varlist = {};

if isempty(varlist0)
    for current_variable_index=1:length(list_of_variables)
        if isequal(list_of_variables(current_variable_index).name,'freq') ...
                || isequal(list_of_variables(current_variable_index).name,'time') ...
                || isequal(list_of_variables(current_variable_index).name,'data') ...
                || isequal(list_of_variables(current_variable_index).name,'varlist') ...
                || isequal(list_of_variables(current_variable_index).name,'varlist0') ...
                || isequal(list_of_variables(current_variable_index).name,'list_of_variables') ...
                || isequal(list_of_variables(current_variable_index).name,'tex') ...
            continue
        end
        if list_of_variables(current_variable_index).global || list_of_variables(current_variable_index).persistent
            % A variable cannot be a global or persistent variable.
            continue
        end
        if list_of_variables(current_variable_index).complex || ~strcmp(list_of_variables(current_variable_index).class,'double')
            % A variable cannot be complex.
            continue
        end
        if list_of_variables(current_variable_index).size(2)>1
            % A variable must be passed as a column vector.
            continue
        end
        try
            eval(['data = [data, ' list_of_variables(current_variable_index).name '];'])
            eval(['varlist = {varlist{:}, ''' list_of_variables(current_variable_index).name '''};'])
        catch
            error(['load_m_file:: All the vectors (variables) in ' inputname(1) ' must have the same number of rows (observations)!'])
        end
    end
else
    for current_variable_index=1:length(varlist0)
       eval(['data = [data, ' varlist0{current_variable_index} '];'])
    end
    varlist = varlist0;
end

%@test:1
%$ % Create a data m-file
%$ fid = fopen('data_m_file.m','w');
%$ fprintf(fid,'FREQ__ = 4;');
%$ fprintf(fid,'INIT__ = ''1938Q4'';');
%$ fprintf(fid,'NAMES__ = {''azert'';''yuiop''};');
%$ fprintf(fid,'TEX__ = {''azert'';''yuiop''};');
%$ fprintf(fid,'azert = [1; 2; 3; 4; 5];');
%$ fprintf(fid,'yuiop = [2; 3; 4; 5; 6];');
%$ fclose(fid);
%$
%$ % Try to read the data m-file
%$ try
%$     datafile = 'data_m_file';
%$     [freq,init,data,varlist,tex] = load_m_file_data(datafile);
%$     t(1) = 1;
%$ catch exception
%$     t(1) = 0;
%$     T = all(t);
%$     LOG = getReport(exception,'extended');
%$     return
%$ end
%$
%$ % Check the results.
%$ t(2) = dassert(freq,4);
%$ t(3) = dassert(isa(init,'dates'),1);
%$ t(4) = dassert(init.freq,4);
%$ t(5) = dassert(init.time,[1938 4]);
%$ t(6) = dassert(varlist,{'azert';'yuiop'});
%$ t(7) = dassert(tex,{'azert';'yuiop'});
%$ t(8) = dassert(data(:,1),[1;2;3;4;5]);
%$ t(9) = dassert(data(:,2),[2;3;4;5;6]);
%$ T = all(t);
%@eof:1
