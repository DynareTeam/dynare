function A = subsasgn(A,S,B) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{A} =} subsasgn (@var{A}, @var{S}, @var{B})
%! @anchor{@dseries/subsasgn}
%! @sp 1
%! Overloads the subsasgn method for the Dynare time series class (@ref{dseries}).
%! @end deftypefn
%@eod:

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

merge_dseries_objects = 1;

switch length(S)
    case 1
      switch S(1).type
        case '{}' % Multiple variable selection.
          if ~isequal(numel(S(1).subs),numel(unique(S(1).subs)))
              error('dseries::subsasgn: Wrong syntax!')
          end
          for i=1:numel(S(1).subs)
              element = S(1).subs{i};
              idArobase = strfind(element,'@');
              if ~isempty(idArobase)
                  switch length(idArobase)
                    case 2
                      idComma = strfind(element(idArobase(1)+1:idArobase(2)-1),',');
                      if ~isempty(idComma)
                          elements = cell(1,numel(idComma)+1); j = 1;
                          expression = element(idArobase(1)+1:idArobase(2)-1);
                          while ~isempty(expression)
                              [token, expression] = strtok(expression,',');
                              elements(j) = {[element(1:idArobase(1)-1), token, element(idArobase(2)+1:end)]};
                              j = j + 1;
                          end
                          S(1).subs = replace_object_in_a_one_dimensional_cell_array(S(1).subs, elements(:), i);
                      else
                          error('dseries::subsasgn: Wrong syntax, matlab''s regular expressions cannot be used here!')
                      end
                    case 4
                      idComma_1 = strfind(element(idArobase(1)+1:idArobase(2)-1),',');
                      idComma_2 = strfind(element(idArobase(3)+1:idArobase(4)-1),',');
                      if ~isempty(idComma_1)
                          elements = cell(1,(numel(idComma_1)+1)*(numel(idComma_2)+1)); j = 1;
                          expression_1 = element(idArobase(1)+1:idArobase(2)-1);
                          while ~isempty(expression_1)
                              [token_1, expression_1] = strtok(expression_1,',');
                              expression_2 = element(idArobase(3)+1:idArobase(4)-1);
                              while ~isempty(expression_2)
                                  [token_2, expression_2] = strtok(expression_2,',');
                                  elements(j) = {[element(1:idArobase(1)-1), token_1, element(idArobase(2)+1:idArobase(3)-1), token_2, element(idArobase(4)+1:end)]};
                                  j = j+1;
                              end
                          end
                          S(1).subs = replace_object_in_a_one_dimensional_cell_array(S(1).subs, elements(:), i);
                      else
                          error('dseries::subsasgn: Wrong syntax, matlab''s regular expressions cannot be used here!')
                      end
                    otherwise
                      error('dseries::subsasgn: Wrong syntax!')
                  end
              end
          end
          if isempty(B)
              for i=1:length(S(1).subs)
                  A = remove(A,S(1).subs{i});
              end
              return
          end
          if ~isequal(length(S(1).subs),vobs(B))
              error('dseries::subsasgn: Wrong syntax!')
          end
          if ~isequal(S(1).subs(:),B.name)
              for i = 1:vobs(B)
                  if ~isequal(S(1).subs{i},B.name{i})
                      % Rename a variable.
                      id = find(strcmp(S(1).subs{i},A.name));
                      if isempty(id)
                          % Add a new variable a change its name.
                          B.name(i) = {S(1).subs{i}};
                          B.tex(i) = {name2tex(S(1).subs{i})};
                      else
                          % Rename variable and change its content.
                          B.name(i) = A.name(id);
                          B.tex(i) = A.tex(id);
                      end
                  end
              end
          end
        case '.'
          if isequal(S(1).subs,'init') && isdates(B) && isequal(length(B),1)
              % Change the initial date (update dates member)
              A.dates = B:B+(nobs(A)-1);
              return
          elseif isequal(S(1).subs,'dates') && isdates(B)
              % Overwrite the dates member
              A.dates = B;
              return
          elseif ismember(S(1).subs,{'data','name','tex'})
              error(['dseries::subsasgn: You cannot overwrite ' S(1).subs ' member!'])
          elseif ~isequal(S(1).subs,B.name)
              % Single variable selection.
              if ~isequal(S(1).subs,B.name{1})
                  % Rename a variable.
                  id = find(strcmp(S(1).subs,A.name));
                  if isempty(id)
                      % Add a new variable a change its name.
                      B.name(1) = {S(1).subs};
                      B.tex(1) = {name2tex(S(1).subs)};
                  else
                      % Rename variable and change its content.
                      B.name(1) = A.name(id);
                      B.tex(1) = A.tex(id);
                  end
              end
          end
        case '()' % Date(s) selection
          if isdates(S(1).subs{1}) || isdate(S(1).subs{1})
                if isdate(S(1).subs{1})
                    Dates = dates(S(1).subs{1});
                else
                    Dates = S(1).subs{1};
                end
              [junk, tdx] = intersect(A.dates.time,Dates.time,'rows');
              if isdseries(B)
                  [junk, tdy] = intersect(B.dates.time,Dates.time,'rows');
                  if isempty(tdy)
                      error('dseries::subsasgn: Periods of the dseries objects on the left and right hand sides must intersect!')
                  end
                  if ~isequal(vobs(A), vobs(B))
                      error('dseries::subsasgn: Dimension error! The number of variables on the left and right hand side must match.')
                  end
                  A.data(tdx,:) = B.data(tdy,:);
                  merge_dseries_objects = 0;
              elseif isnumeric(B)
                  merge_dseries_objects = 0;
                  if isequal(length(tdx),rows(B))
                      if isequal(columns(A.data),columns(B))
                          A.data(tdx,:) = B;
                      else
                          error('dseries::subsasgn: Dimension error! The number of variables on the left and right hand side must match.')
                      end
                  else
                      error('dseries::subsassgn: Dimension error! The number of periods on the left and right hand side must match.')
                  end
              else
                  error('dseries::subsasgn: The object on the right hand side must be a dseries object or a numeric array!')
              end
          elseif ischar(S(1).subs{1}) && isequal(S(1).subs{1},':') && isempty(A)
              if isnumeric(B)
                  if isequal(rows(B),1)
                      A.data = repmat(B,A.dates.ndat,1);
                  elseif isequal(rows(B),A.dates.ndat)
                      A.data = B;
                  else
                      error('dseries::subsasgn: Wrong syntax!')
                  end
                  if isempty(A.name)
                      A.name = default_name(vobs(A));
                      A.tex = name2tex(A.name);
                  end
              elseif isdseries(B)
                  if isequal(nobs(B), 1)
                      A.data = repmat(B.data,A.dates.ndat,1);
                  elseif isequal(nobs(B), A.dates.ndat)
                      A.data = B;
                  else
                      error('dseries::subsasgn: Wrong syntax!')
                  end
                  if isempty(A.name)
                      A.name = B.name;
                      A.tex = B.tex;
                  end
              end
              return
          else
              error('dseries::subsasgn: Wrong syntax!')
          end
        otherwise
          error('dseries::subsasgn: Wrong syntax!')
      end
  case 2
    merge_dseries_objects = 0;
    if ((isequal(S(1).type,'{}') || isequal(S(1).type,'.')) && isequal(S(2).type,'()'))
        if isequal(S(1).type,'{}')
            sA = extract(A,S(1).subs{:});
        else
            sA = extract(A,S(1).subs);
        end
        if (isdseries(B) && isequal(vobs(sA), vobs(B))) || (isnumeric(B) && isequal(vobs(sA),columns(B))) || (isnumeric(B) && isequal(columns(B),1))
            if isdates(S(2).subs{1})
                [junk, tdx] = intersect(sA.dates.time,S(2).subs{1}.time,'rows');
                if isdseries(B)
                    [junk, tdy] = intersect(B.dates.time,S(2).subs{1}.time,'rows');
                    if isempty(tdy)
                        error('dseries::subsasgn: Periods of the dseries objects on the left and right hand sides must intersect!')
                    end
                    sA.data(tdx,:) = B.data(tdy,:);
                elseif isnumeric(B)
                    merge_dseries_objects = 0;
                    if isequal(length(tdx),rows(B))
                        if isequal(columns(sA.data),columns(B))
                            sA.data(tdx,:) = B;
                        elseif isequal(size(B,2),1)
                            sA.data(tdx,:) = repmat(B,1,columns(sA.data));
                        else
                            error('dseries::subsasgn: Dimension error! The number of variables on the left and right hand side must match.')
                        end
                    else
                        if isequal(columns(sA.data),columns(B)) && isequal(rows(B),1)
                            sA.data(tdx,:) = repmat(B,length(tdx),1);
                        elseif isequal(rows(B),1)
                            sA.data(tdx,:) = B;
                        else
                            error('dseries::subsassgn: Dimension error! The number of periods on the left and right hand side must match.')
                        end
                    end
                else
                    error('dseries::subsasgn: The object on the right hand side must be a dseries object or a numeric array!')
                end
            else
                error('dseries::subsasgn: Wrong syntax!')
            end
            A = merge(A,sA);
        else
            error('dseries::subsasgn: Dimension error! The number of variables on the left and right hand side must match.')
        end
    end
  otherwise
    error('dseries::subsasgn: Wrong syntax!')
end

if isempty(A)
    % Assign variables to an empty dseries object.
    A = B;
    return
end

if merge_dseries_objects
    A = merge(A,B);
end

%@test:1
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A2'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[A(:,1), B, A(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$
%$ % Apply the exponential function to the second variable.
%$ ts1{'A2'} = ts1{'A2'}.exp;
%$
%$ % Instantiate a time series object.
%$
%$    t(1) = dassert(ts1.vobs,3);
%$    t(2) = dassert(ts1.nobs,10);
%$    t(3) = dassert(ts1.name{2},'A2');
%$    t(4) = dassert(ts1.name{1},'A1');
%$    t(5) = dassert(ts1.name{3},'A3');
%$    t(6) = dassert(ts1.data,[A(:,1), exp(A(:,2)), A(:,3)],1e-15);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables.
%$ ts1{'A1'} = ts1{'A1'}.log;
%$ ts1{'A3'} = ts1{'A3'}.log;
%$
%$ % Instantiate a time series object.
%$
%$    t(1) = dassert(ts1.vobs,3);
%$    t(2) = dassert(ts1.nobs,10);
%$    t(3) = dassert(ts1.name{2},'A2');
%$    t(4) = dassert(ts1.name{1},'A1');
%$    t(5) = dassert(ts1.name{3},'A3');
%$    t(6) = dassert(ts1.data,[log(A(:,1)), A(:,2), log(A(:,3))],1e-15);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables of ts1.
%$ ts1{'A1','A3'} = ts1{'A1','A3'}.log;
%$
%$ t(1) = dassert(ts1.vobs,3);
%$ t(2) = dassert(ts1.nobs,10);
%$ t(3) = dassert(ts1.name{1},'A1') && dassert(ts1.name{2},'A2') && dassert(ts1.name{3},'A3');
%$ t(4) = dassert(ts1.data,[log(A(:,1)), A(:,2), log(A(:,3))],1e-15);
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'A1';'B2';'B3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables.
%$ ts1.A1 = ts2.A1;
%$
%$ % Instantiate a time series object.
%$
%$ t(1) = dassert(ts1.vobs,3);
%$ t(2) = dassert(ts1.nobs,10);
%$ t(3) = dassert(ts1.name{1},'A1');
%$ t(4) = dassert(ts1.data(:,1),B(:,1), 1e-15);
%$ t(5) = dassert(ts1.data(:,2:3),A(:,2:3), 1e-15);
%$ T = all(t);
%@eof:5

%@test:6
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,2);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'B1';'B2'},[]);
%$
%$ % Call tested routine.
%$ try
%$     ts1.B2 = ts2.B2.log;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,4);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name{1},'A1');
%$    t(5) = dassert(ts1.name{2},'A2');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.name{4},'B2');
%$    t(8) = dassert(ts1.data,[A(:,1), A(:,2), A(:,3), log(B(:,2))],1e-15);
%$ end
%$ T = all(t);
%@eof:6

%@test:7
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,2);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'B1';'B2'},[]);
%$
%$ % Append B2 to the first object.
%$ ts1{'B2'} = ts2{'B2'};
%$ t(1) = 1;
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,4);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name{1},'A1');
%$    t(5) = dassert(ts1.name{2},'A2');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(6) = dassert(ts1.name{4},'B2');
%$    t(7) = dassert(ts1.data,[A(:,1), A(:,2), A(:,3), B(:,2)],1e-15);
%$ end
%$ T = all(t);
%@eof:7

%@test:8
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A4'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,4);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.name{4},'A4');
%$    t(8) = dassert(ts1.data,[A, B],1e-15);
%$ end
%$ T = all(t);
%@eof:8

%@test:9
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,2);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'A1';'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A1','A4'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,4);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.name{4},'A4');
%$    t(8) = dassert(ts1.data,[B(:,1), A(:,2:3), B(:,2)],1e-15);
%$ end
%$ T = all(t);
%@eof:9


%@test:10
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,[],{'A1';'B1';'B2'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A@1,2@','A4'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,4);
%$    t(3) = dassert(ts1.nobs,10);
%$    t(4) = dassert(ts1.name{1},'A1');
%$    t(5) = dassert(ts1.name{2},'A2');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.name{4},'A4');
%$    t(8) = dassert(ts1.data,[B(:,1:2), A(:,3), B(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:10

%@test:11
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,5);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,[],{'A_1';'A_2';'A_3'},[]);
%$ ts2 = dseries(B,[],{'A_1';'A_2';'B_1';'B_2';'B_3'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'@A,B@_@1,2@'} = ts2{'@A,B@_@1,2@'};
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    %t(2) = dassert(ts1.vobs,4);
%$    %t(3) = dassert(ts1.nobs,10);
%$    %t(4) = dassert(ts1.name,{'A1','A2';'A3';'B1';'B2'});
%$    %t(5) = dassert(ts1.data,[B(:,1:2), A(:,3), B(:,3:4)],1e-15);
%$ end
%$ T = all(t);
%@eof:11

%@test:12
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1{'A1'}(rg) = ts2{'B1'}(rg);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); B(3:7); A(8:end,1)], A(:,2:3)],1e-15);
%$ end
%$ T = all(t);
%@eof:12

%@test:13
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1{'A1'}(rg) = B(3:7);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); B(3:7); A(8:end,1)], A(:,2:3)],1e-15);
%$ end
%$ T = all(t);
%@eof:13

%@test:14
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1.A1(rg) = B(3:7);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); B(3:7); A(8:end,1)], A(:,2:3)],1e-15);
%$ end
%$ T = all(t);
%@eof:14

%@test:15
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1.A1(rg) = sqrt(pi);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); repmat(sqrt(pi),5,1); A(8:end,1)], A(:,2:3)],1e-15);
%$ end
%$ T = all(t);
%@eof:15

%@test:16
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1{'A1','A2'}(rg) = sqrt(pi);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); repmat(sqrt(pi),5,1); A(8:end,1)], [A(1:2,2); repmat(sqrt(pi),5,1); A(8:end,2)], A(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:16

%@test:17
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1{'A1','A2'}(rg) = [sqrt(pi), pi];
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); repmat(sqrt(pi),5,1); A(8:end,1)], [A(1:2,2); repmat(pi,5,1); A(8:end,2)], A(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:17

%@test:18
%$ % Define a datasets.
%$ A = rand(40,3); B = rand(40,1);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$ ts2 = dseries(B,'1950Q1',{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     d1 = dates('1950Q3');
%$     d2 = dates('1951Q3');
%$     rg = d1:d2;
%$     ts1{'A1','A2'}(rg) = ones(5,1);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,[[A(1:2,1); ones(5,1); A(8:end,1)], [A(1:2,2); ones(5,1); A(8:end,2)], A(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:18

%@test:19
%$ % Define a datasets.
%$ A = rand(40,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$
%$ % Instantiate a dates object.
%$ dd = dates('1952Q1');
%$
%$ % modify first object.
%$ try
%$     ts1.init = dd;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,A,1e-15);
%$    t(8) = dassert(ts1.init,dd);
%$    t(9) = dassert(ts1.dates(1),dd);
%$ end
%$ T = all(t);
%@eof:19

%@test:20
%$ % Define a datasets.
%$ A = rand(40,3);
%$
%$ % Instantiate two dseries object.
%$ ts1 = dseries(A,'1950Q1',{'A1';'A2';'A3'},[]);
%$
%$ % Instantiate a dates object.
%$ dd = dates('1952Q1');
%$
%$ % modify first object.
%$ try
%$     ts1.dates = dd:(dd+(ts1.nobs-1));
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts1.vobs,3);
%$    t(3) = dassert(ts1.nobs,40);
%$    t(4) = dassert(ts1.name{2},'A2');
%$    t(5) = dassert(ts1.name{1},'A1');
%$    t(6) = dassert(ts1.name{3},'A3');
%$    t(7) = dassert(ts1.data,A,1e-15);
%$    t(8) = dassert(ts1.init,dd);
%$    t(9) = dassert(ts1.dates(1),dd);
%$ end
%$ T = all(t);
%@eof:20

%@test:21
%$ % Define a datasets.
%$ A = rand(4,3);
%$
%$ % Instantiate an empty dseries object.
%$ ts = dseries(dates('1950Q1'):dates('1950Q4'));
%$
%$ % Populate ts
%$ try
%$     ts(:) = A;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts.vobs,3);
%$    t(3) = dassert(ts.nobs,4);
%$    t(4) = dassert(ts.data,A,1e-15);
%$ end
%$ T = all(t);
%@eof:21

%@test:22
%$ % Define a datasets.
%$ A = rand(1,3);
%$
%$ % Instantiate an empty dseries object.
%$ ts = dseries(dates('1950Q1'):dates('1950Q4'));
%$
%$ % Populate ts
%$ try
%$     ts(:) = A;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dassert(ts.vobs,3);
%$    t(3) = dassert(ts.nobs,4);
%$    t(4) = dassert(ts.data,repmat(A,4,1),1e-15);
%$ end
%$ T = all(t);
%@eof:22

%@test:23
%$ % Instantiate a dseries object.
%$ ts0 = dseries(randn(10,6), '1999y');
%$
%$ % Try to remove Variable_2 and Variable_3
%$ try
%$     ts0{'Variable_@2,3@'} = [];
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$    % Try to display Variable_2 and Variable_3 again (should fail because already removed)
%$    try
%$       ts0{'Variable_@2,3@'};
%$       t(2) = 0;
%$    catch
%$       t(2) = 1;
%$    end
%$ end
%$ T = all(t);
%@eof:23