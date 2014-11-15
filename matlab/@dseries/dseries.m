function ts = dseries(varargin) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{ts} =} dseries (@var{a},@var{b},@var{c},@var{d})
%! @anchor{dseries}
%! @sp 1
%! Constructor for the Dynare time series class.
%! @sp 2
%! @strong{Inputs}
%! @sp 2
%! If @code{nargin==0} then an empty dseries object is created. The object can be populated with data subsequently using the overloaded subsref method.
%! @sp 2
%! If @code{nargin==1} and if the input argument is a @ref{dates} object, then a dseries object without data is created. This object can be populated with the overload subsref method.
%! @sp 2
%! If @code{nargin==1} and if the input argument is a string for the name of a csv, m or mat file containing data, then a dseries object is created from these data.
%! @sp 2
%! If @code{nargin>1}:
%! @sp 1
%! @table @ @var
%! @item a
%! T*1 vector or T*N matrix of data.
%! @item b
%! Initial date. For Quaterly, Monthly or Weekly data, b must be a string. For yearly data or if the frequence is not defined b must be an integer.
%! @item c
%! N*1 cell array of strings or N*q array of strings. Names of the N time series.
%! @item d
%! N*1 cell array of strings or N*p array of characters. TeX names of the N time series.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Dynare time series object.
%! @end table
%! @sp 2
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item data
%! Array of doubles (nobs*vobs).
%! @item nobs
%! Scalar integer, the number of observations.
%! @item vobs
%! Scalar integer, the number of variables.
%! @item name
%! Cell array of strings, names of the variables.
%! @item tex
%! Cell array of strings, tex names of the variables.
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item init
%! @ref{dates} object, initial date of the dataset.
%! @end table
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

if nargin>0 && ischar(varargin{1}) && isequal(varargin{1},'initialize')
    ts = struct;
    ts.data  = [];
    ts.name  = {};
    ts.tex   = {};
    ts.dates = dates();
    ts = class(ts,'dseries');
    assignin('base','emptydseriesobject',ts);
    return
end

ts = evalin('base','emptydseriesobject');

switch nargin
  case 0
    %  Create an empty dseries object.
    return
  case 1
    if isdates(varargin{1})
        switch length(varargin{1})
          case 0
            error(['dseries::dseries: Input (identified as a dates object) must be non empty!'])
          case 1
            % Create an empty dseries object with an initial date.
            ts.dates = varargin{1};
          otherwise
            % A range of dates is passed to the constructor
            ts.dates = varargin{1};
        end
        return
    elseif ischar(varargin{1})
        % Create a dseries object loading data in a file (*.csv, *.m, *.mat).
        if isempty(varargin{1})
            error('dseries:: Wrong calling sequence! Input argument cannot be an empty string.')
        elseif check_file_extension(varargin{1},'m')
            [freq,init,data,varlist,tex] = load_m_file_data(varargin{1});
        elseif check_file_extension(varargin{1},'mat')
            [freq,init,data,varlist,tex] = load_mat_file_data(varargin{1});
        elseif check_file_extension(varargin{1},'csv')
            [freq,init,data,varlist] = load_csv_file_data(varargin{1});
            tex = [];
        elseif check_file_extension(varargin{1},'xls') || check_file_extension(varargin{1},'xlsx')
            if isglobalinbase('options_')
                % Check that the object is instantiated within a dynare session so that options_ global structure exists.
                % Should provide latter a mechanism to pass range and sheet to dseries constructor...
                range = evalin('base','options_.xls_range');
                sheet = evalin('base','options_.xls_sheet');
            else
                % By default only the (whole) first sheet is loaded.
                range = [];
                sheet = [];
            end
            [freq,init,data,varlist] = load_xls_file_data(varargin{1}, sheet, range);
            tex = [];
        else
            error(['dseries:: I''m not able to load data from ' varargin{1} '!'])
        end
        ts.data = data;
        ts.name = varlist;
        ts.dates = init:init+(nobs(ts)-1);
        if isempty(tex)
            ts.tex = name2tex(varlist);
        else
            ts.tex = tex;
        end
    elseif isnumeric(varargin{1}) && isequal(ndims(varargin{1}),2)
        ts.data = varargin{1};
        ts.name = default_name(vobs(ts));
        ts.tex = name2tex(ts.name);
        ts.dates = dates(1,1):dates(1,1)+(nobs(ts)-1);
    end
  case {2,3,4}
    if isequal(nargin,2) && ischar(varargin{1}) && isdates(varargin{2})
        % Instantiate dseries object with a data file and force the initial date to be as given by the second input argument.
        ds = dseries(varargin{1});
        ts = dseries(ds.data, varargin{2}, ds.name, ds.tex);
        return
    end
    if isequal(nargin,2) && ischar(varargin{1}) && ischar(varargin{2}) && isdate(varargin{2})
        % Instantiate dseries object with a data file and force the initial date to be as given by the second input argument.
        ds = dseries(varargin{1});
        ts = dseries(ds.data, dates(varargin{2}), ds.name, ds.tex);
        return
    end
    a = varargin{1};
    b = varargin{2};
    if nargin<4
        d = {};
    else
        d = varargin{4};
        if ~iscell(d) && ~isempty(d)
            d = cellstr(d);
        end
    end
    if nargin<3
        c = {};
    else
        c = varargin{3};
        if ~iscell(c) && ~isempty(c)
            c = cellstr(c);
        end
    end
    % Get data, number of observations and number of variables.
    ts.data = a;
    % Get the first date and set the frequency.
    if isempty(b)
        init = dates(1,1);
    elseif (isdates(b) && isequal(length(b),1))
        init = b;
    elseif ischar(b) && isdate(b)% Weekly, Monthly, Quaterly or Annual data (string).
        init = dates(b);
    elseif (isnumeric(b) && isscalar(b) && isint(b)) % Yearly data.
        init = dates([num2str(b) 'Y']);
    elseif isdates(b) % Range of dates
        init = b(1);
        if nobs(ts)>1 && ~isequal(b.ndat,nobs(ts))
            message =   'dseries::dseries: If second input is a range, its number of elements must match ';
            message = char(message, '                  the number of rows in the first input, unless the first input');
            message = char(message, '                  has only one row.');
            skipline()
            disp(message);
            error(' ');
        elseif isequal(nobs(ts), 1)
            ts.data = repmat(ts.data,b.ndat,1);
        end
        ts.dates = b;
    elseif (isnumeric(b) && isint(b)) % Range of yearly dates.
        message = 'dseries::dseries: Not implemented! If you need to define a range of years';
        message = char(message, '                  you have to pass a dates object as the second input argument.');
        disp(message)
        error(' ')
    else
        error('dseries::dseries: Wrong calling sequence!');
    end
    % Get the names of the variables.
    if ~isempty(c)
        if vobs(ts)==length(c)
            for i=1:vobs(ts)
                ts.name = vertcat(ts.name, c(i));
            end
        else
            error('dseries::dseries: The number of declared names does not match the number of variables!')
        end
    else
        ts.name = default_name(vobs(ts));
    end
    if ~isempty(d)
        if vobs(ts)==length(d)
            for i=1:vobs(ts)
                ts.tex = vertcat(ts.tex, d(i));
            end
        else
            error('dseries::dseries: The number of declared tex names does not match the number of variables!')
        end
    else
        ts.tex = name2tex(ts.name);
    end
  otherwise
    error('dseries::dseries: Can''t instantiate the class, wrong calling sequence!')
end

if isempty(ts.dates)
    ts.dates = init:init+(nobs(ts)-1);
end

%@test:1
%$ % Test if we can instantiate an empty dseries object.
%$ try
%$     ts = dseries();
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(4,1);
%$
%$ try
%$     aa = dates('1938M11');
%$     ts = dseries(aa);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,12);
%$     t(3) = dassert(ts.init.freq,12);
%$     t(4) = dassert(ts.init.time,[1938, 11]);
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dseries('dynseries_test_data.m');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1994, 3]);
%$     t(5) = dassert(ts.vobs,2);
%$     t(6) = dassert(ts.nobs,100);
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dseries('dynseries_test_data.mat');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1994, 3]);
%$     t(5) = dassert(ts.vobs,2);
%$     t(6) = dassert(ts.nobs,100);
%$ end
%$
%$ T = all(t);
%@eof:4

%@test:5
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dseries('dynseries_test_data.csv');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1990, 1]);
%$     t(5) = dassert(ts.vobs,4);
%$     t(6) = dassert(ts.nobs,4);
%$     t(7) = dassert(ts.name,{'azert';'yuiop';'qsdfg';'jklm'});
%$     t(8) = dassert(ts.tex,{'azert';'yuiop';'qsdfg';'jklm'});
%$ end
%$
%$ T = all(t);
%@eof:5

%@test:6
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dseries(transpose(1:5),[]);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,1);
%$     t(3) = dassert(ts.init.freq,1);
%$     t(4) = dassert(ts.init.time,[1, 1]);
%$     t(5) = dassert(ts.vobs,1);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Variable_1'});
%$     t(8) = dassert(ts.tex,{'Variable\\_1'});
%$ end
%$
%$ T = all(t);
%@eof:6

%@test:7
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dseries(transpose(1:5),'1950Q1');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1950, 1]);
%$     t(5) = dassert(ts.vobs,1);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Variable_1'});
%$     t(8) = dassert(ts.tex,{'Variable\\_1'});
%$ end
%$
%$ T = all(t);
%@eof:7

%@test:8
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dseries([transpose(1:5), transpose(6:10)],'1950q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1950, 1]);
%$     t(5) = dassert(ts.vobs,2);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Output'; 'Consumption'});
%$     t(8) = dassert(ts.tex,{'Y_t'; 'C_t'});
%$ end
%$
%$ T = all(t);
%@eof:8

%@test:9
%$ try
%$     ts = dseries('dynseries_test_data-1.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1990, 1]);
%$     t(5) = dassert(ts.vobs,3);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'GDP';'Consumption';'CPI'});
%$     t(8) = dassert(ts.tex,{'GDP';'Consumption';'CPI'});
%$ end
%$
%$ T = all(t);
%@eof:9

%@test:10
%$ try
%$     ts = dseries('dynseries_test_data-2.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1990, 1]);
%$     t(5) = dassert(ts.vobs,3);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Variable_1';'Variable_2';'Variable_3'});
%$     t(8) = dassert(ts.tex,{'Variable\\_1';'Variable\\_2';'Variable\\_3'});
%$ end
%$
%$ T = all(t);
%@eof:10

%@test:11
%$ try
%$     ts = dseries('dynseries_test_data-3.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(ts.freq,1);
%$     t(3) = dassert(ts.init.freq,1);
%$     t(4) = dassert(ts.init.time,[1, 1]);
%$     t(5) = dassert(ts.vobs,3);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'Variable_1';'Variable_2';'Variable_3'});
%$     t(8) = dassert(ts.tex,{'Variable\\_1';'Variable\\_2';'Variable\\_3'});
%$ end
%$
%$ T = all(t);
%@eof:11

%@test:12
%$ try
%$     ts = dseries('dynseries_test_data-4.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(ts.freq,1);
%$     t(3) = dassert(ts.init.freq,1);
%$     t(4) = dassert(ts.init.time,[1, 1]);
%$     t(5) = dassert(ts.vobs,3);
%$     t(6) = dassert(ts.nobs,5);
%$     t(7) = dassert(ts.name,{'GDP';'Consumption';'CPI'});
%$     t(8) = dassert(ts.tex,{'GDP';'Consumption';'CPI'});
%$ end
%$
%$ T = all(t);
%@eof:12

%@test:13
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dseries(transpose(1:4),dates('1990Q1'):dates('1990Q4'));
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1990, 1]);
%$     t(5) = dassert(ts.vobs,1);
%$     t(6) = dassert(ts.nobs,4);
%$ end
%$
%$ T = all(t);
%@eof:13

%@test:14
%$ t = zeros(7,1);
%$
%$ try
%$     ts = dseries([1, 2],dates('1990Q1'):dates('1990Q4'));
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dassert(ts.freq,4);
%$     t(3) = dassert(ts.init.freq,4);
%$     t(4) = dassert(ts.init.time,[1990, 1]);
%$     t(5) = dassert(ts.vobs,2);
%$     t(6) = dassert(ts.nobs,4);
%$     t(7) = dassert(ts.data, [ones(4,1), 2*ones(4,1)]);
%$ end
%$
%$ T = all(t);
%@eof:14

%@test:15
%$ try
%$     ts = dseries([1; 2],dates('1990Q1'):dates('1990Q4'));
%$     t = 0;
%$ catch
%$     t = 1;
%$ end
%$
%$ T = all(t);
%@eof:15
