function B = subsref(A, S) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{us} =} subsref (@var{ts},S)
%! @anchor{@dseries/subsref}
%! @sp 1
%! Overloads the subsref method for the Dynare time series class (@ref{dseries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Dynare time series object instantiated by @ref{dseries}.
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
%! or a dseries object built by extracting a variable from @var{ts}, or a dseries object containing a subsample of the all the variable in @var{ts}.
%! @end table
%! @sp 2
%! @strong{Example 1.} Let @var{ts} be a dseries object containing three variables named 'A1', 'A2' and 'A3'. Then the following syntax:
%! @example
%!   us = ts.A1;
%! @end example
%!will create a new dseries object @var{us} containing the variable 'A1'.
%! @sp 1
%! @strong{Example 2.} Let @var{ts} be a dseries object. Then the following syntax:
%! @example
%!   us = ts.log;
%! @end example
%!will create a new dseries object @var{us} containing all the variables of @var{ts} transformed by the neperian logarithm.
%! @sp 1
%! @strong{Example 3.} Let @var{ts} be a dseries object. The following syntax:
%! @example
%!   us = ts(3:50);
%! @end example
%!will create a new dseries object @var{us} by selecting a subsample out of @var{ts}.
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2014 Dynare Team
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
      case {'data','name','tex','dates'}        % Public members.
        if length(S)>1 && isequal(S(2).type,'()') && isempty(S(2).subs)
            error(['dseries::subsref: ' S(1).subs ' is not a method but a member!'])
        end
        B = builtin('subsref', A, S(1));
      case {'log','exp','ygrowth','qgrowth','ydiff','qdiff','abs'}         % Give "dot access" to public methods without args.
        B = feval(S(1).subs,A);
        if length(S)>1 && isequal(S(2).type,'()') && isempty(S(2).subs)
            S = shiftS(S,1);
        end
      case 'nobs'
        % Returns the number of observations.
        B = rows(A.data);
      case 'vobs'
        % Returns the number of variables.
        B = columns(A.data);
      case 'init'
        % Returns a dates object (first date).
        B = A.dates(1);
      case 'last'
        % Returns a dates object (last date).
        B = A.dates(end);
      case 'freq'
        % Returns an integer characterizing the data frequency (1, 4, 12 or 52)
        B = A.dates.freq;
      case {'lag','lead','hptrend','hpcycle','chain','detrend','exist'} % Methods with less than two arguments.
        if length(S)>1 && isequal(S(2).type,'()')
            if isempty(S(2).subs)
                B = feval(S(1).subs,A);
                S = shiftS(S,1);
            else
                if ~ischar(S(2).subs{1}) && length(S(2).subs{1})>1
                    error(['dseries::subsref: ' S(1).subs{1} ' method admits no more than one argument!'])
                end
                B = feval(S(1).subs,A,S(2).subs{1});
                S = shiftS(S,1);
            end
        else
            B = feval(S(1).subs,A);
        end
      case {'cumsum','insert','pop','cumprod','remove'} % Methods with less than three argument.
        if length(S)>1 && isequal(S(2).type,'()')
            if isempty(S(2).subs)
                B = feval(S(1).subs,A);
                S = shiftS(S,1);
            else
                if length(S(2).subs)>2
                    error(['dseries::subsref: ' S(1).subs{1} ' method admits no more than two arguments!'])
                end
                B = feval(S(1).subs,A,S(2).subs{:});
                S = shiftS(S,1);
            end
        else
            B = feval(S(1).subs,A);
        end
      case 'baxter_king_filter'
        if length(S)>1 && isequal(S(2).type,'()')
            if isempty(S(2).subs)
                B = feval(S(1).subs,A);
                S = shiftS(S,1);
            else
                B = feval(S(1).subs,A,S(2).subs{1})
                S = shiftS(S,1);
            end
        else
            B = feval(S(1).subs,A);
        end
      case 'save'                                                        % Save dseries object on disk (default is a csv file).
        B = NaN;
        if isequal(length(S),2)
            if strcmp(S(2).type,'()')
                if isempty(S(2).subs)
                    save(A,inputname(1));
                else
                    if isempty(S(2).subs{1})
                        save(A,inputname(1),S(2).subs{2});
                    else
                        save(A,S(2).subs{:});
                    end
                end
                S = shiftS(S,1);
            else
                error('dseries::subsref: Wrong syntax.')
            end
        elseif isequal(length(S),1)
            save(A,inputname(1));
        else
            error('dseries::subsref: Call to save method must come in last position!')
        end
      case 'size'
        if isequal(length(S),2) && strcmp(S(2).type,'()')
            if isempty(S(2).subs)
                [x,y] = size(A);
                B = [x, y];
            else
                B = size(A,S(2).subs{1});
            end
            S = shiftS(S,1);
        elseif isequal(length(S),1)
            [x,y] = size(A);
            B = [x, y];
        else
            error('dseries::subsref: Call to size method must come in last position!')
        end
      case {'set_names','rename','tex_rename'}
        B = feval(S(1).subs,A,S(2).subs{:});
        S = shiftS(S,1);
      case {'disp'}
        feval(S(1).subs,A);
        return
      otherwise                                                            % Extract a sub-object by selecting one variable.
        ndx = find(strcmp(S(1).subs,A.name));
        if ~isempty(ndx)
            B = dseries();
            B.data = A.data(:,ndx);
            B.name = A.name(ndx);
            B.tex = A.tex(ndx);
            B.tex  = deblank(A.tex(ndx,:));
            B.dates = A.dates;
        else
            error('dseries::subsref: Unknown public method, public member or variable!')
        end
    end
  case '()'
    if ischar(S(1).subs{1}) && ~isdate(S(1).subs{1})
        % If ts is an empty dseries object, populate this object by reading data in a file.
        if isempty(A)
            B = dseries(S(1).subs{1});
        else
            error(['dseries::subsref: dseries object ''' inputname(1) '''  is not empty!'])
        end
    elseif isa(S(1).subs{1},'dynTimeIndex')
        % shift backward/forward (lag/lead) dseries object
        shift = S(1).subs{1}.index;
        if shift>0
            B = feval('lead',A,shift);
        elseif shift<0
            B = feval('lag',A,-shift);
        else
            % Do nothing.
            B = A;
        end
    elseif isscalar(S(1).subs{1}) && isnumeric(S(1).subs{1}) && isint(S(1).subs{1})
        % Input is also interpreted as a backward/forward operator
        if S(1).subs{1}>0
            B = feval('lead', A, S(1).subs{1});
        elseif S(1).subs{1}<0
            B = feval('lag', A, -S(1).subs{1});
        else
            % Do nothing.
            B = A;
        end
    elseif isdates(S(1).subs{1}) || isdate(S(1).subs{1})
        if isdate(S(1).subs{1})
            Dates = dates(S(1).subs{1});
        else
            Dates = S(1).subs{1};
        end
        % Test if Dates is out of bounds
        if min(Dates)<min(A.dates)
            error(['dseries::subsref: Indices are out of bounds! Subsample cannot start before ' date2string(A.dates(1)) '.'])
        end
        if  max(Dates)>max(A.dates)
            error(['dseries::subsref: Indices are out of bounds! Subsample cannot end after ' date2string(A.dates(end)) '.'])
        end
        % Extract a subsample using a dates object
        [junk,tdx] = intersect(A.dates.time,Dates.time,'rows');
        B = dseries();
        B.data = A.data(tdx,:);
        B.name = A.name;
        B.tex  = A.tex;
        B.dates = A.dates(tdx);
    elseif isvector(S(1).subs{1}) && all(isint(S(1).subs{1}))
        error('dseries::subsref: It is not possible to select observations with a vector of integers. You have to index with a dates object instead!');
    else
        error('dseries::subsref: I have no idea of what you are trying to do!')
    end
  case '{}'
    if iscellofchar(S(1).subs)
        B = extract(A,S(1).subs{:});
    elseif isequal(length(S(1).subs),1) && all(isint(S(1).subs{1}))
        idx = S(1).subs{1};
        if max(idx)>size(A.data,2) || min(idx)<1
            error('dseries::subsref: Indices are out of bounds!')
        end
        B = dseries();
        B.data = A.data(:,idx);
        B.name = A.name(idx);
        B.tex  = A.tex(idx);
        B.dates = A.dates;
    else
        error('dseries::subsref: What the Hell are you tryin'' to do?!')
    end
  otherwise
    error('dseries::subsref: What the Hell are you doin'' here?!')
end

S = shiftS(S,1);
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
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1(ts1.dates(2:9));
%$
%$ % Expected results.
%$ e.data = [transpose(2:9),2*transpose(2:9)];
%$ e.nobs = 8;
%$ e.vobs = 2;
%$ e.name = {'A1';'A2'};
%$ e.freq = 1;
%$ e.init = dates(1,2);
%$
%$ % Check the results.
%$ t(1) = dassert(a.data,e.data);
%$ t(2) = dassert(a.nobs,e.nobs);
%$ t(3) = dassert(a.vobs,e.vobs);
%$ t(4) = dassert(a.freq,e.freq);
%$ t(5) = dassert(a.init,e.init);
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
%$ ts1 = dseries(A,[],A_name,[]);
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
%$ e.init = dates(1,1);
%$
%$ % Check the results.
%$ t(1) = dassert(a.data,e.data);
%$ t(2) = dassert(a.init,e.init);
%$ t(3) = dassert(a.nobs,e.nobs);
%$ t(4) = dassert(a.vobs,e.vobs);
%$ t(5) = dassert(a.freq,e.freq);
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
%$ ts1 = dseries(A,[],A_name,[]);
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
%$ e.init = dates(1,1);
%$
%$ % Check the results.
%$ t(1) = dassert(a.data,e.data);
%$ t(2) = dassert(a.nobs,e.nobs);
%$ t(3) = dassert(a.vobs,e.vobs);
%$ t(4) = dassert(a.freq,e.freq);
%$ t(5) = dassert(a.init,e.init);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Create an empty dseries object.
%$ dataset = dseries();
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
%$     t(2) = dassert(A.nobs,4);
%$     t(3) = dassert(A.vobs,4);
%$     t(4) = dassert(A.freq,4);
%$     t(5) = dassert(A.init,dates('1990Q1'));
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
%$ ts1 = dseries(A,[],A_name,[]);
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
%$ e.init = dates(1,1);
%$
%$ t(1) = dassert(e.data,a.data);
%$ t(2) = dassert(e.nobs,a.nobs);
%$ t(3) = dassert(e.vobs,a.vobs);
%$ t(4) = dassert(e.name,a.name);
%$ t(5) = dassert(e.init,a.init);
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
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1{'GDP_[0-9]'};
%$ b = ts1{'[A-Z]_1$'};
%$
%$ % Expected results.
%$ e1.data = A(:,1:12);
%$ e1.nobs = 10;
%$ e1.vobs = 12;
%$ e1.name = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; 'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'};
%$ e1.freq = 1;
%$ e1.init = dates(1,1);
%$ e2.data = A(:,[1 13]);
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
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts1.save('ts1');
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
%$    ts1 = dseries(A,[],A_name,[]);
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
%$ ts1 = dseries(A,'1971Q1',A_name,[]);
%$
%$ % Define the range of a subsample.
%$ range = dates('1971Q2'):dates('1971Q4');
%$ % Call the tested method.
%$ a = ts1(range);
%$
%$ % Expected results.
%$ e.data = A(2:4,:);
%$ e.nobs = 3;
%$ e.vobs = 3;
%$ e.name = {'A1';'A2';'B1'};
%$ e.freq = 4;
%$ e.init = dates('1971Q2');
%$
%$ t(1) = dassert(e.data,a.data);
%$ t(2) = dassert(e.nobs,a.nobs);
%$ t(3) = dassert(e.vobs,a.vobs);
%$ t(4) = dassert(e.name,a.name);
%$ t(5) = dassert(e.init,a.init);
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
%$ ts1 = dseries(A,'1971Q1',A_name,[]);
%$
%$ % Test the size method.
%$ B = ts1.size();
%$ C = ts1.size(1);
%$ D = ts1.size(2);
%$ E = ts1.size;
%$
%$ t(1) = dassert(B,[60, 3]);
%$ t(2) = dassert(E,[60, 3]);
%$ t(3) = dassert(C,60);
%$ t(4) = dassert(D,3);
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
%$ ts1 = dseries(A,'1971Q1',A_name,[]);
%$
%$ % Test the size method.
%$ B = ts1{1};
%$ C = ts1{[1,3]};
%$ D = ts1{'A1'};
%$
%$ t(1) = dassert(B.name{1},'A1');
%$ t(2) = dassert(B.data,A(:,1));
%$ t(3) = dassert(C.name{1},'A1');
%$ t(4) = dassert(C.data(:,1),A(:,1));
%$ t(5) = dassert(C.name{2},'B1');
%$ t(6) = dassert(C.data(:,2),A(:,3));
%$ t(7) = dassert(D.name{1},'A1');
%$ t(8) = dassert(D.data,A(:,1));
%$ T = all(t);
%@eof:11

%@test:12
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    if isoctave
%$        ts1.save('ts1');
%$    else
%$        ts1.save();
%$    end
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:12

%@test:13
%$ try
%$     data = transpose(0:1:50);
%$     ts = dseries(data,'1950Q1');
%$     a = ts.lag;
%$     b = ts.lead;
%$     tt = dynTimeIndex();
%$     c = ts(tt-1);
%$     d = ts(tt+1);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)>1
%$     t(2) = (a==c);
%$     t(3) = (b==d);
%$ end
%$
%$ T = all(t);
%@eof:13

%@test:14
%$ try
%$     data = transpose(0:1:50);
%$     ts = dseries(data,'1950Q1');
%$     a = ts.lag;
%$     b = ts.lead;
%$     c = ts(-1);
%$     d = ts(1);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)>1
%$     t(2) = (a==c);
%$     t(3) = (b==d);
%$ end
%$
%$ T = all(t);
%@eof:14

%@test:15
%$ try
%$     ds = dseries(transpose(1:5));
%$     ts = ds(ds.dates(2:3));
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)>1
%$     t(2) = isdseries(ts);
%$     t(3) = isequal(ts.data,ds.data(2:3));
%$ end
%$
%$ T = all(t);
%@eof:15

%@test:16
%$ try
%$     ds = dseries(transpose(1:5));
%$     ts = ds(ds.dates(2:6));
%$     t(1) = 0;
%$ catch
%$     t(1) = 1;
%$ end
%$
%$ T = all(t);
%@eof:16

%@test:17
%$ try
%$     ds = dseries(transpose(1:5));
%$     ts = ds(dates('1Y'):dates('6Y'));
%$     t(1) = 0;
%$ catch
%$     t(1) = 1;
%$ end
%$
%$ T = all(t);
%@eof:17

%@test:18
%$ try
%$     ds = dseries(transpose(1:5));
%$     ts = ds(dates('-2Y'):dates('4Y'));
%$     t(1) = 0;
%$ catch
%$     t(1) = 1;
%$ end
%$
%$ T = all(t);
%@eof:18