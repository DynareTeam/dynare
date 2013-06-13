function B = subsref(A, S)
%@info:
%! @deftypefn {Function File} {@var{us} =} subsref (@var{ts},S)
%! @anchor{@dynSeries/subsref}
%! @sp 1
%! Overloads the subsref method for the Dynare time series class (@ref{dynSeries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Dynare time series object instantiated by @ref{dynSeries}.
%! @item S
%! Matlab's structure array S with two fields, type and subs. The type field is string containing '()', '@{@}', or '.', where '()' specifies
%! integer subscripts, '@{@}' specifies cell array subscripts, and '.' specifies subscripted structure fields. The subs field is a cell array
%! or a string containing the actual subscripts (see matlab's documentation).
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item us
%! Dynare time series object. Depending on the calling sequence @var{us} is a transformation of @var{ts} obtained by applying a public method on @var{ts},
%! or a dynSeries object built by extracting a variable from @var{ts}, or a dynSeries object containing a subsample of the all the variable in @var{ts}.
%! @end table
%! @sp 2
%! @strong{Example 1.} Let @var{ts} be a dynSeries object containing three variables named 'A1', 'A2' and 'A3'. Then the following syntax:
%! @example
%!   us = ts.A1;
%! @end example
%!will create a new dynSeries object @var{us} containing the variable 'A1'.
%! @sp 1
%! @strong{Example 2.} Let @var{ts} be a dynSeries object. Then the following syntax:
%! @example
%!   us = ts.log;
%! @end example
%!will create a new dynSeries object @var{us} containing all the variables of @var{ts} transformed by the neperian logarithm.
%! @sp 1
%! @strong{Example 3.} Let @var{ts} be a dynSeries object. The following syntax:
%! @example
%!   us = ts(3:50);
%! @end example
%!will create a new dynSeries object @var{us} by selecting a subsample out of @var{ts}.
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2013 Dynare Team
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

switch S(1).type
  case '.'
    switch S(1).subs
      case {'data','nobs','vobs','name','tex','freq','time','init'}        % Public members.
        B = builtin('subsref', A, S(1));
      case {'log','exp','ygrowth','qgrowth','ydiff','qdiff','lag'}         % Give "dot access" to public methods.
        B = feval(S(1).subs,A);
      case {'save'}                                                        % Save dynSeries object on disk (default is a csv file). 
        B = NaN;
        if length(S)==2 && strcmp(S(2).type,'()')
            save(A,S(2).subs{:});
            S = shiftS(S);
        else
            save(A);
        end
      case {'size'}
        if length(S)==2 && strcmp(S(2).type,'()') && ~isempty(S(2).subs)
            B = size(A,S(2).subs{1});
            S = shiftS(S);
        else
            [x,y] = size(A);
            B = [x, y];
        end
      case {'rename','tex_rename'}
        B = feval(S(1).subs,A,S(2).subs{:});
        S = shiftS(S);
      otherwise                                                            % Extract a sub-object by selecting one variable.
        ndx = strmatch(S(1).subs,A.name,'exact');
        if ~isempty(ndx)
            B = dynSeries();
            B.data = A.data(:,ndx);
            B.name = A.name(ndx);
            B.tex = A.tex(ndx);
            B.tex  = deblank(A.tex(ndx,:));
            B.nobs = A.nobs;
            B.vobs = 1;
            B.freq = A.freq;
            B.init = A.init;
            B.time = A.time;
        else
            error('dynSeries::subsref: Unknown public method, public member or variable!')
        end
    end    
  case '()'
    if ischar(S(1).subs{1})
        % If ts is an empty dynSeries object, populate this object by reading data in a file.
        if isempty(A)
            B = dynSeries(S(1).subs{1});
        else
            error(['dynSeries::subsref: dynSeries object ''' inputname(1) '''  is not empty!'])
        end
    elseif isa(S(1).subs{1},'dynDates')
        % Extract a subsample using a dynDates object
        [junk,tdx] = intersect(A.time.time,S(1).subs{1}.time,'rows');
        B = dynSeries();
        B.data = A.data(tdx,:);
        B.name = A.name;
        B.tex  = A.tex;
        B.nobs = length(tdx);
        B.vobs = A.vobs;
        B.freq = A.freq;
        B.init = A.init+tdx(1);
        B.time = A.time(tdx,:);
    elseif isvector(S(1).subs{1}) && all(isint(S(1).subs{1}))
        % Extract a subsample using a vector of integers (observation index).
        if all(S(1).subs{1}>0) && all(S(1).subs{1}<=A.nobs)
            if size(A.data,2)>1
                S(1).subs = [S(1).subs, ':'];
            end
            B.data = builtin('subsref', A.data, S(1));
            B.nobs = size(B.data,1);
            B.vobs = A.vobs;
            B.freq = A.freq;
            B.time = builtin('subsref', A.time, S(1));
            B.init = A.init+S(1).subs{1}(1);
            B.name = A.name;
            B.tex  = A.tex;
        else
            error('dynSeries::subsref: Indices are out of bounds!')
        end
    else
        error('dynSeries::subsref: I have no idea of what you are trying to do!')
    end
  case '{}'
    if ischar(S(1).subs{1})
        S(1).subs = strtrim(mat2cell(S(1).subs{1}, ones(1, size(S(1).subs{1},1)), size(S(1).subs{1},2)));
    end
    if iscellofchar(S(1).subs)
        B = extract(A,S(1).subs{:});
    elseif isequal(length(S(1).subs),1) && all(isint(S(1).subs{1}))
        idx = S(1).subs{1};
        if max(idx)>A.vobs || min(idx)<1
            error('dynSeries::subsref: Indices are out of bounds!')
        end
        B = dynSeries();
        B.data = A.data(:,idx);
        B.name = A.name(idx);
        B.tex  = A.tex(idx);
        B.nobs = A.nobs;
        B.vobs = length(idx);
        B.freq = A.freq;
        B.init = A.init;
        B.time = A.time;
    else
        error('dynSeries::subsref: What the Hell are you tryin'' to do?!')
    end
  otherwise
    error('dynSeries::subsref: What the Hell are you doin'' here?!')
end

S = shiftS(S);
if ~isempty(S)
    B = subsref(B, S);
end

%@test:1
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1(2:9);
%$
%$ % Expected results.
%$ e.data = [transpose(2:9),2*transpose(2:9)];
%$ e.nobs = 8;
%$ e.vobs = 2;
%$ e.name = {'A1';'A2'};
%$ e.freq = 1;
%$ e.init = dynDate(2);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.nobs,e.nobs);
%$ t(3) = dyn_assert(a.vobs,e.vobs);
%$ t(4) = dyn_assert(a.freq,e.freq);
%$ t(5) = dyn_assert(a.init,e.init);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1.A1;
%$
%$ % Expected results.
%$ e.data = transpose(1:10);
%$ e.nobs = 10;
%$ e.vobs = 1;
%$ e.name = {'A1'};
%$ e.freq = 1;
%$ e.init = dynDate(1);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.init,e.init);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1.log;
%$
%$ % Expected results.
%$ e.data = log(A);
%$ e.nobs = 10;
%$ e.vobs = 2;
%$ e.name = {'A1';'A2'};
%$ e.freq = 1;
%$ e.init = dynDate(1);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.nobs,e.nobs);
%$ t(3) = dyn_assert(a.vobs,e.vobs);
%$ t(4) = dyn_assert(a.freq,e.freq);
%$ t(5) = dyn_assert(a.init,e.init);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Create an empty dynSeries object.
%$ dataset = dynSeries();
%$
%$ t = zeros(5,1);
%$
%$ try
%$    A = dataset('dynseries_test_data.csv');
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ % Check the results.
%$ if length(t)>1
%$     t(2) = dyn_assert(A.nobs,4);
%$     t(3) = dyn_assert(A.vobs,4);
%$     t(4) = dyn_assert(A.freq,4);
%$     t(5) = dyn_assert(A.init,dynDate('1990Q1'));
%$ end
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10),3*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2';'B1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1{'A1','B1'};
%$
%$ % Expected results.
%$ e.data = A(:,[1,3]);
%$ e.nobs = 10;
%$ e.vobs = 2;
%$ e.name = {'A1';'B1'};
%$ e.freq = 1;
%$ e.init = dynDate(1);
%$
%$ t(1) = dyn_assert(e.data,a.data);
%$ t(2) = dyn_assert(e.nobs,a.nobs);
%$ t(3) = dyn_assert(e.vobs,a.vobs);
%$ t(4) = dyn_assert(e.name,a.name);
%$ t(5) = dyn_assert(e.init,a.init);
%$ T = all(t);
%@eof:5

%@test:6
%$ % Define a data set.
%$ A = rand(10,24);
%$
%$ % Define names
%$ A_name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; 'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'; 'HICP_1';'HICP_2';'HICP_3'; 'HICP_4'; 'HICP_5'; 'HICP_6'; 'HICP_7'; 'HICP_8'; 'HICP_9'; 'HICP_10'; 'HICP_11'; 'HICP_12';};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1{'GDP_@0-9@'};
%$ b = ts1{'@A-Z@_1'};
%$
%$ % Expected results.
%$ e1.data = A(:,1:12);
%$ e1.nobs = 10;
%$ e1.vobs = 12;
%$ e1.name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; 'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'};
%$ e1.freq = 1;
%$ e1.init = dynDate(1);
%$ e2.data = A(:,[1, 13]);
%$ e2.nobs = 10;
%$ e2.vobs = 2;
%$ e2.name = {'GDP_1';'HICP_1'};
%$ e2.freq = 1;
%$ e2.init = dynDate(1);
%$
%$ % Check results.
%$ t(1) = dyn_assert(e1.data,a.data);
%$ t(2) = dyn_assert(e1.nobs,a.nobs);
%$ t(3) = dyn_assert(e1.vobs,a.vobs);
%$ t(4) = dyn_assert(e1.name,a.name);
%$ t(5) = dyn_assert(e1.init,a.init);
%$ t(6) = dyn_assert(e2.data,b.data);
%$ t(7) = dyn_assert(e2.nobs,b.nobs);
%$ t(8) = dyn_assert(e2.vobs,b.vobs);
%$ t(9) = dyn_assert(e2.name,b.name);
%$ t(10) = dyn_assert(e2.init,b.init);
%$ T = all(t);
%@eof:6

%@test:7
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    ts1.save;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:7

%@test:8
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    ts1.save('test_generated_data_file','m');
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:8

%@test:9
%$ % Define a data set.
%$ A = [transpose(1:60),2*transpose(1:60),3*transpose(1:60)];
%$
%$ % Define names
%$ A_name = {'A1';'A2';'B1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,'1971Q1',A_name,[]);
%$
%$ % Define the range of a subsample.
%$ range = dynDate('1971Q2'):dynDate('1971Q4');
%$ % Call the tested method.
%$ a = ts1(range);
%$
%$ % Expected results.
%$ e.data = A(2:4,:);
%$ e.nobs = 3;
%$ e.vobs = 3;
%$ e.name = {'A1';'A2';'B1'};
%$ e.freq = 4;
%$ e.init = dynDate('1971Q2');
%$
%$ t(1) = dyn_assert(e.data,a.data);
%$ t(2) = dyn_assert(e.nobs,a.nobs);
%$ t(3) = dyn_assert(e.vobs,a.vobs);
%$ t(4) = dyn_assert(e.name,a.name);
%$ t(5) = dyn_assert(e.init,a.init);
%$ T = all(t);
%@eof:9

%@test:10
%$ % Define a data set.
%$ A = [transpose(1:60),2*transpose(1:60),3*transpose(1:60)];
%$
%$ % Define names
%$ A_name = {'A1';'A2';'B1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,'1971Q1',A_name,[]);
%$
%$ % Test the size method.
%$ B = ts1.size();
%$ C = ts1.size(1);
%$ D = ts1.size(2);
%$ E = ts1.size;
%$
%$ t(1) = dyn_assert(B,[60, 3]);
%$ t(2) = dyn_assert(E,[60, 3]);
%$ t(3) = dyn_assert(C,60);
%$ t(4) = dyn_assert(D,3);
%$ T = all(t);
%@eof:10

%@test:11
%$ % Define a data set.
%$ A = [transpose(1:60),2*transpose(1:60),3*transpose(1:60)];
%$
%$ % Define names
%$ A_name = {'A1';'A2';'B1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,'1971Q1',A_name,[]);
%$
%$ % Test the size method.
%$ B = ts1{1};
%$ C = ts1{[1,3]};
%$ D = ts1{'A1'};
%$
%$ t(1) = dyn_assert(B.name{1},'A1');
%$ t(2) = dyn_assert(B.data,A(:,1));
%$ t(3) = dyn_assert(C.name{1},'A1');
%$ t(4) = dyn_assert(C.data(:,1),A(:,1));
%$ t(5) = dyn_assert(C.name{2},'B1');
%$ t(6) = dyn_assert(C.data(:,2),A(:,3));
%$ t(7) = dyn_assert(D.name{1},'A1');
%$ t(8) = dyn_assert(D.data,A(:,1));
%$ T = all(t);
%@eof:11
