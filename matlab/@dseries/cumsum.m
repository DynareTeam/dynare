function B = cumsum(varargin) % --*-- Unitary tests --*--

% Overloads matlab's cumsum function for dseries objects.
%
% INPUTS 
%  o A     dseries object [mandatory].
%  o d     dates object [optional]
%  o v     dseries object with one observation [optional]
%
% OUTPUTS 
%  o B     dseries object.

% Copyright (C) 2013 Dynare Team
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

% Get indices of the columns without NaNs
idx = find(~any(isnan(varargin{1}.data)));

if isempty(idx)
    error('dseries::cumsum: All the variables have NaNs. The cumulated sum cannot be computed!')
end

if ~isequal(idx(:),transpose(1:vobs(varargin{1})))
    warning('dseries::cumsum: The cumulated sum is not computed for some variables because they have NaNs!')
end

switch nargin
    case 1
      % Initialize the output.
      B = varargin{1};
      % Perform the cumulated sum
      B.data(:,idx) = cumsum(B.data(:,idx));
      % Change the name of the variables
      for i=1:vobs(B)
          B.name(i) = {['cumsum(' B.name{i} ')']};
          B.tex(i) = {['\sum_t ' B.tex{i}]};
      end
  case 2
    if isdseries(varargin{2})
        if ~isequal(vobs(varargin{1}), vobs(varargin{2}))
            error('dseries::cumsum: First and second input arguments must be dseries objects with the same number of variables!')
        end
        if ~isequal(varargin{1}.name, varargin{2}.name)
            warning('dseries::cumsum: First and second input arguments must be dseries objects do not have the same variables!')
        end
        if ~isequal(nobs(varargin{2}),1)
            error('dseries::cumsum: Second input argument must be a dseries object with only one observation!')
        end
        B = cumsum(varargin{1});
        B.data = bsxfun(@plus,B.data,varargin{2}.data);
    elseif isdates(varargin{2})
        B = cumsum(varargin{1});
        t = find(B.dates==varargin{2});
        if isempty(t)
            if varargin{2}==(firstdate(B)-1)
                return
            else
                error(['dseries::cumsum: date ' date2string(varargin{2}) ' is not in the sample!'])
            end
        end
        B.data = bsxfun(@minus,B.data,B.data(t,:));
    else
        error('dseries::cumsum: Second input argument must be a dseries object or a dates object!')
    end
  case 3
    if ~isdates(varargin{2})
        error('dseries::cumsum: Second input argument must be a dates object!')
    end
    if ~isdseries(varargin{3})
        error('dseries::cumsum: Third input argument must be a dseries object!')
    end
    if ~isequal(vobs(varargin{1}), vobs(varargin{3}))
        error('dseries::cumsum: First and third input arguments must be dseries objects with the same number of variables!')
    end
    if ~isequal(varargin{1}.name, varargin{3}.name)
        warning('dseries::cumsum: First and third input arguments must be dseries objects do not have the same variables!')
    end
    if ~isequal(nobs(varargin{3}),1)
        error('dseries::cumsum: Third input argument must be a dseries object with only one observation!')
    end
    B = cumsum(varargin{1});
    t = find(B.dates==varargin{2});
    if isempty(t)
        if varargin{2}==(firstdate(B)-1)
            B.data = bsxfun(@plus,B.data,varargin{3}.data);
            return
        else
            error(['dseries::cumsum: date ' date2string(varargin{2}) ' is not in the sample!'])
        end
    end
    B.data = bsxfun(@plus,B.data,varargin{3}.data-B.data(t,:));
  otherwise
    error('dseries::cumsum: Wrong number of input arguments!')
end

%@test:1
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ ts1 = cumsum(ts1);
%$
%$ % Expected results.
%$ ts2 = dseries(transpose(1:10), [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dassert(ts1,ts2);
%$ warning_config
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ ts1 = ts1.cumsum();
%$
%$ % Expected results.
%$ ts2 = dseries(transpose(1:10), [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dassert(ts1,ts2);
%$ warning_config
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$ 
%$ % Call the tested method.
%$ ts1 = cumsum(ts1,dates('3Y'));
%$
%$ % Expected results.
%$ ts2 = dseries([-2; -1; 0; 1; 2; 3; 4; 5; 6; 7], [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dassert(ts1,ts2);
%$ warning_config
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$ ts2 = dseries(pi, [], A_name, []);
%$
%$ % Call the tested method.
%$ ts3 = cumsum(ts1,dates('3Y'),ts2);
%$
%$ % Expected results.
%$ ts4 = dseries([-2; -1; 0; 1; 2; 3; 4; 5; 6; 7]+pi, [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dassert(ts3,ts4);
%$ warning_config
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dseries(A,[],A_name,[]);
%$ ts2 = dseries(pi, [], A_name, []);
%$
%$ % Call the tested method.
%$ ts3 = ts1.cumsum(dates('3Y'),ts2);
%$
%$ % Expected results.
%$ ts4 = dseries([-2; -1; 0; 1; 2; 3; 4; 5; 6; 7]+pi, [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dassert(ts3,ts4);
%$ warning_config
%$ T = all(t);
%@eof:5