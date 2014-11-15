function B = subsref(A,S) % --*-- Unitary tests --*--

% Overload the subsref method for dates objects.
%
% INPUTS 
%  o A     dates object.
%  o S     matlab's structure.
%
% OUTPUTS 
%  o B     dates object.
%
% REMARKS 
%  See the matlab's documentation about the subsref method.

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
      case {'time','freq','ndat'}% Access public members.
        if length(S)>1 && isequal(S(2).type,'()') && isempty(S(2).subs)
            error(['dates::subsref: ' S(1).subs ' is not a method but a member!'])
        end
        B = builtin('subsref', A, S(1));
      case {'sort','unique','double','isempty','length','char'}% Public methods (without input arguments)
        B = feval(S(1).subs,A);
        if length(S)>1 && isequal(S(2).type,'()') && isempty(S(2).subs)
           S = shiftS(S,1);
        end
      case {'append','pop'}% Public methods (with arguments).
        if isequal(S(2).type,'()')
            B = feval(S(1).subs,A,S(2).subs{:});
            S = shiftS(S,1);
        else
            error('dates::subsref: Something is wrong in your syntax!')
        end
      case {'disp'}
        feval(S(1).subs,A);
        return
      otherwise
        error('dates::subsref: Unknown public member or method!')
    end
  case '()'
    if isempty(A)
        if isempty(A.freq)
            % Populate an empty dates object with time member (freq is not specified). Needs two or three inputs. First input is an integer
            % scalar specifying the frequency. Second input is either the time member (a n*2 array of integers) or a column vector with n
            % elements (the first column of the time member --> years). If the the second input is a row vector and if A.freq~=1 a third input
            % is necessary. The third input is n*1 vector of integers between 1 and A.freq (second column of the time member --> subperiods).
            B = dates();
            % First input is the frequency.
            if isfreq(S(1).subs{1})
                if ischar(S(1).subs{1})
                    B.freq = string2freq(S(1).subs{1});
                else
                    B.freq = S(1).subs{1};
                end
            else
                error('dates::subsref: First input must be a frequency!')
            end
            if isequal(length(S(1).subs),2)
                % If two inputs are provided, the second input must be a n*2 array except if frequency is annual.
                [n, m] = size(S(1).subs{2});
                if m>2
                    error('dates::subsref: Second input argument array cannot have more than two columns!')
                end
                if ~isequal(m,2) && ~isequal(B.freq,1)
                    error('dates::subsref: Second argument has to be a n*2 array for non annual frequency!')
                end
                if ~all(all(S(1).subs{2}))
                    error('dates::subsref: Second argument has be an array of intergers!')
                end
                if m>1 && ~issubperiod(S(1).subs{2}(:,2),B.freq)
                    error(['dates::subsref: Elements in the second column of the first input argument are not legal subperiods (should be integers betwwen 1 and ' num2str(B.freq) ')!'])
                end
                if isequal(m,2)
                    B.time = S(1).subs{2};
                elseif isequal(m,1)
                    B.time = [S(1).subs{2}, ones(n,1)];
                else
                    error('dates::subsref: This is a bug!')
                end
                B.ndat = rows(B.time);
            elseif isequal(length(S(1).subs),3)
                % If three inputs are provided, the second and third inputs are column verctors of integers (years and subperiods).
                if ~iscolumn(S(1).subs{2}) && ~all(isint(S(1).subs{2}))
                    error('dates::subsref: Second input argument must be a column vector of integers!')
                end
                n1 = size(S(1).subs{2},1);
                if ~iscolumn(S(1).subs{3}) && ~issubperiod(S(1).subs{3}, B.freq)
                    error(['dates::subsref: Third input argument must be a column vector of subperiods (integers between 1 and ' num2str(B.freq) ')!'])
                end
                n2 = size(S(1).subs{3},1);
                if ~isequal(n1,n2)
                    error('dates::subsref: Second and third input arguments must have the same number of elements!')
                end
                B.time = [S(1).subs{2}, S(1).subs{3}];
                B.ndat = rows(B.time);
            else
                error('dates::subsref: Wrong calling sequence!')
            end
        else
            % Populate an empty dates object with time member (freq is already specified).
            % Needs one (time) or two (first and second columns of time for years and subperiods) inputs.
            B = A;
            if isequal(length(S(1).subs),2)
                if ~iscolumn(S(1).subs{1}) && ~all(isint(S(1).subs{1}))
                    error('dates::subsref: First argument has to be a column vector of integers!')
                end
                n1 = size(S(1).subs{1},1);
                if ~iscolumn(S(1).subs{2}) && ~issubperiod(S(1).subs{2}, B.freq)
                    error(['dates::subsref: Second argument has to be a column vector of subperiods (integers between 1 and ' num2str(B.freq) ')!'])
                end
                n2 = size(S(1).subs{2},1);
                if ~isequal(n2,n1)
                    error('dates::subsref: First and second argument must have the same number of rows!')
                end
                B.time = [S(1).subs{1}, S(1).subs{2}];
                B.ndat = rows(B.time);
            elseif isequal(length(S(1).subs),1)
                [n, m] = size(S(1).subs{1});
                if ~isequal(m,2) && ~isequal(B.freq,1)
                    error('dates::subsref: First argument has to be a n*2 array!')
                end
                if ~all(isint(S(1).subs{1}(:,1)))
                    error('dates::subsref: First column of the first argument has to be a column vector of integers!')
                end
                if m>1 && issubperiod(S(1).subs{1}(:,1), B.freq)
                    error(['dates::subsref: The second column of the first input argument has to be a column  vector of subperiods (integers between 1 and ' num2str(B.freq) ')!'])
                end
                if isequal(m,2)
                    B.time = S(1).subs{1};
                elseif isequal(m,1) && isequal(B.freq,1)
                    B.time = [S(1).subs{1}, ones(n,1)];
                else
                    error('dates::subsref: This is a bug!')
                end
                B.ndat = rows(B.time);
            else
                error('dates::subsref: Wrong number of inputs!')
            end
        end
    else
        % dates object A is not empty. We extract some dates
        if isvector(S(1).subs{1}) && all(isint(S(1).subs{1})) && all(S(1).subs{1}>0) && all(S(1).subs{1}<=A.ndat)
            B = dates();
            B.freq = A.freq;
            B.time = A.time(S(1).subs{1},:);
            B.ndat = rows(B.time);
        else
            error(['dates::subsref: indices has to be a vector of positive integers less than or equal to ' int2str(A.ndat) '!'])
        end
    end
  otherwise
    error('dates::subsref: Something is wrong in your syntax!')
end

S = shiftS(S,1);
if ~isempty(S)
    B = subsref(B, S);
end

%@test:1
%$ % Define a dates object
%$ B = dates('1950Q1','1950Q2','1950Q3','1950Q4','1951Q1');
%$
%$ % Try to extract a sub-dates object.
%$ d = B(2:3);
%$
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(d.freq,B.freq);
%$     t(3) = dassert(d.time,[1950 2; 1950 3]);
%$     t(4) = dassert(d.ndat,2);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a dates object
%$ B = dates('1950Q1'):dates('1960Q3');
%$
%$ % Try to extract a sub-dates object and apply a method
%$ 
%$ d = B(2:3).sort ;
%$
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(d.freq,B.freq);
%$     t(3) = dassert(d.time,[1950 2; 1950 3]);
%$     t(4) = dassert(d.ndat,2);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a dates object
%$ B = dates('1950Q1'):dates('1960Q3');
%$
%$ % Try to extract a sub-dates object and apply a method
%$
%$ d = B(2:3).sort() ;
%$
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(d.freq,B.freq);
%$     t(3) = dassert(d.time,[1950 2; 1950 3]);
%$     t(4) = dassert(d.ndat,2);
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a dates object
%$ B = dates('1950Q1','1950Q2','1950Q3','1950Q4','1951Q1');
%$
%$ % Try to extract a sub-dates object.
%$ d = B(2);
%$
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dassert(d.freq,B.freq);
%$     t(3) = dassert(d.time,[1950 2]);
%$     t(4) = dassert(d.ndat,1);
%$ end
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define an empty dates object with quaterly frequency.
%$ qq = dates('Q');
%$
%$ % Define a ranges of dates using qq.
%$ try
%$     r1 = qq(1950,1):qq(1950,3);%qq([1950, 3]);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$ if t(1)
%$     try
%$         r2 = qq([1950, 1; 1950, 2; 1950, 3]);
%$         t(2) = 1;
%$     catch
%$         t(2) = 0;
%$     end
%$ end
%$ if t(1) && t(2)
%$     try
%$         r3 = qq(1950*ones(3,1), transpose(1:3));
%$         t(3) = 1;
%$     catch
%$         t(3) = 0;
%$     end
%$ end
%$
%$ if t(1) && t(2) && t(3)
%$     t(4) = dassert(r1,r2);
%$     t(5) = dassert(r1,r3);
%$ end
%$ T = all(t);
%@eof:5

%@test:6
%$ % Define an empty dates object with quaterly frequency.
%$ date = dates();
%$
%$ % Define a ranges of dates using qq.
%$ try
%$     r1 = date(4,1950,1):date(4,[1950, 3]);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$ if t(1)
%$     try
%$         r2 = date(4,[1950, 1; 1950, 2; 1950, 3]);
%$         t(2) = 1;
%$     catch
%$         t(2) = 0;
%$     end
%$ end
%$ if t(1) && t(2)
%$     try
%$         r3 = date(4,1950*ones(3,1), transpose(1:3));
%$         t(3) = 1;
%$     catch
%$         t(3) = 0;
%$     end
%$ end
%$
%$ if t(1) && t(2) && t(3)
%$     t(4) = dassert(r1,r2);
%$     t(5) = dassert(r1,r3);
%$ end
%$ T = all(t);
%@eof:6

