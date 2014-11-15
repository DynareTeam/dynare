function A = extract(B,varargin) % --*-- Unitary tests --*--

% Extract some variables from a database.
    
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

A = dseries();

% Get the names of the variables to be extracted from dseries object B.
VariableName_ = {};
for i=1:nargin-1
    VariableName = varargin{i};
    idArobase = strfind(VariableName,'@');
    if mod(length(idArobase),2)
        error('dseries::extract: (Implicit loops) The number of @ symbols must be even!')
    end
    idBracket.open = strfind(VariableName,'[');
    idBracket.close = strfind(VariableName,']');
    if ~isequal(length(idBracket.open),length(idBracket.open))
        error('dseries::extract: (Matlab/Octave''s regular expressions) Check opening and closing square brackets!')
    end
    if length(idArobase)
        NumberOfImplicitLoops = .5*length(idArobase);
        idComma = cell(NumberOfImplicitLoops,1);
        expressions = cell(NumberOfImplicitLoops,1);
        for i=0:NumberOfImplicitLoops-1
            idComma(i+1) = { strfind(VariableName(idArobase(2*i+1)+1:idArobase(2*i+2)-1),',') };
            expressions(i+1) = { VariableName(idArobase(2*i+1)+1:idArobase(2*i+2)-1) };
        end
        if any(cellfun(@isempty,idComma))
            error('dseries::extract: (Implicit loops) Wrong syntax!')
        end
        switch NumberOfImplicitLoops
          case 1
            expression = expressions{1};
            idVariables_ = [];
            while ~isempty(expression)
                [token, expression] = strtok(expression,',');             
                candidate = [VariableName(1:idArobase(1)-1), token, VariableName(idArobase(2)+1:end)];
                id = find(strcmp(candidate,B.name));
                if isempty(id)
                    error(['dseries::extract: (Implicit loops) Variable ''' candidate ''' does not exist in dseries object ''' inputname(1) '''!'])
                else
                    idVariables_ = [idVariables_; id];
                end
            end
            VariableName = B.name(idVariables_);
          case 2
            idVariables_ = [];
            expression_1 = expressions{1};
            while ~isempty(expression_1)
                [token_1, expression_1] = strtok(expression_1,',');
                expression_2 = expressions{2};
                while ~isempty(expression_2)
                    [token_2, expression_2] = strtok(expression_2,',');
                    candidate = [VariableName(1:idArobase(1)-1), token_1, VariableName(idArobase(2)+1:idArobase(3)-1),  token_2, VariableName(idArobase(4)+1:end)];
                    id = find(strcmp(candidate,B.name));
                    if isempty(id)
                        error(['dseries::extract: (Implicit loops) Variable ''' candidate ''' does not exist in dseries object ''' inputname(1) '''!'])
                    else
                        idVariables_ = [idVariables_; id];
                    end
                end
            end
            VariableName = B.name(idVariables_);
          otherwise
            error('dseries::extract: (Implicit loops) Cannot unroll more than two implicit loops!')
        end
        VariableName_ = vertcat(VariableName_,VariableName);
    elseif length(idBracket.open)
        % Matlab/Octave's regular expressions.
        first_block_id = 0;
        last_block_id = 0;
        idVariables = find(isnotempty_cell(regexp(B.name,VariableName,'match')));
        if isempty(idVariables)
            error(['dseries::extract: Can''t find any variable matching ' VariableName ' pattern!'])
        end
        idVariables_ = [];
        for j = 1:length(idVariables)
            first_block_flag = 0;
            if (first_block_id && strcmp(B.name{idVariables(j)}(1:first_block_id),VariableName(1:first_block_id))) || ~first_block_id
                first_block_flag = 1;
            end
            last_block_flag = 0;
            if (last_block_id && strcmp(B.name{idVariables(j)}(end-last_block_id:end),VariableName(end-last_block_id:end))) || ~last_block_id
                last_block_flag = 1;
            end
            if first_block_flag && last_block_flag
                idVariables_ = [idVariables_; idVariables(j)];
            end
        end
        VariableName = B.name(idVariables_);
        VariableName_ = vertcat(VariableName_,VariableName);
    else
        VariableName_ = varargin(:);
    end
end

% Remove trailing white spaces if any
VariableName_ = strtrim(VariableName_);

% Get indices of the selected variables
idVariableName = NaN(length(VariableName_),1);
for i = 1:length(idVariableName)
    idx = find(strcmp(VariableName_{i},B.name));
    if isempty(idx)
        error(['dseries::extract: Variable ' VariableName_{i} ' is not a member of ' inputname(1) '!'])
    end
    idVariableName(i) = idx;
end

A.data = B.data(:,idVariableName);
A.dates = B.dates;
A.name = B.name(idVariableName);
A.tex = B.tex(idVariableName);

function b = isnotempty_cell(CellArray)
    CellArrayDimension = size(CellArray);
    b = NaN(CellArrayDimension);
    for i=1:CellArrayDimension(1)
        for j = 1:CellArrayDimension(2)
            b(i,j) = ~isempty(CellArray{i,j});
        end
    end
    
    
%@test:1
%$ % Define a data set.
%$ A = rand(10,24);
%$
%$ % Define names
%$ A_name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; 'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'; 'HICP_1';'HICP_2';'HICP_3'; 'HICP_4'; 'HICP_5'; 'HICP_6'; 'HICP_7'; 'HICP_8'; 'HICP_9'; 'HICP_10'; 'HICP_11'; 'HICP_12';};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1{'GDP_@1,2,3,4,5@'};
%$ b = ts1{'@GDP,HICP@_1'};
%$
%$ % Expected results.
%$ e1.data = A(:,1:5);
%$ e1.nobs = 10;
%$ e1.vobs = 5;
%$ e1.name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'};
%$ e1.freq = 1;
%$ e1.init = dates(1,1);
%$ e2.data = A(:,[1, 13]);
%$ e2.nobs = 10;
%$ e2.vobs = 2;
%$ e2.name = {'GDP_1';'HICP_1'};
%$ e2.freq = 1;
%$ e2.init = dates(1,1);
%$
%$ % Check results.
%$ t(1) = dassert(e1.data,a.data);
%$ t(2) = dassert(e1.nobs,a.nobs);
%$ t(3) = dassert(e1.vobs,a.vobs);
%$ t(4) = dassert(e1.name,a.name);
%$ t(5) = dassert(e1.init,a.init);
%$ t(6) = dassert(e2.data,b.data);
%$ t(7) = dassert(e2.nobs,b.nobs);
%$ t(8) = dassert(e2.vobs,b.vobs);
%$ t(9) = dassert(e2.name,b.name);
%$ t(10) = dassert(e2.init,b.init);
%$ T = all(t);
%@eof:1


%@test:2
%$ % Define a data set.
%$ A = rand(10,24);
%$
%$ % Define names
%$ A_name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; 'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'; 'HICP_1';'HICP_2';'HICP_3'; 'HICP_4'; 'HICP_5'; 'HICP_6'; 'HICP_7'; 'HICP_8'; 'HICP_9'; 'HICP_10'; 'HICP_11'; 'HICP_12';};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ try
%$   a = ts1{'GDP_@1,2,3,4,55@'};
%$   t = 0;
%$ catch
%$   t = 1;
%$ end
%$
%$ T = all(t);
%@eof:2


%@test:3
%$ % Define a data set.
%$ A = rand(10,24);
%$
%$ % Define names
%$ A_name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; 'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'; 'HICP_1';'HICP_2';'HICP_3'; 'HICP_4'; 'HICP_5'; 'HICP_6'; 'HICP_7'; 'HICP_8'; 'HICP_9'; 'HICP_10'; 'HICP_11'; 'HICP_12';};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ try
%$   a = ts1{'@GDP,HICP@_@1,2,3,4,5@'};
%$   t = 1;
%$ catch
%$   t = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dassert(a.name,{'GDP_1';'GDP_2';'GDP_3';'GDP_4';'GDP_5';'HICP_1';'HICP_2';'HICP_3';'HICP_4';'HICP_5'});
%$ end
%$
%$ T = all(t);
%@eof:3
