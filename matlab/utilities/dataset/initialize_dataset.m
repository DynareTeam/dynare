function dataset_ = initialize_dataset(datafile,varobs,first,nobs,transformation,prefilter,xls)
% Initializes a structure describing the dataset.

% Copyright (C) 2011 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

if isempty(datafile)
    error('Estimation:: You have to declare a dataset file!')
end

if isempty(varobs)
    error('Estimation:: You have to declare a set of observed variables')
end

% Get raw data.
rawdata = read_variables(datafile,varobs,[],xls.sheet,xls.range);

% Get the (default) number of observations.
if isempty(nobs)
    nobs = rows(rawdata)-first+1;
end

% Get the (default) prefilter option.
if isempty(prefilter)
    prefilter = 0;
end

% Fill the dataset structure
dataset_.info.ntobs = nobs;
dataset_.info.nvobs = rows(varobs);
dataset_.info.varobs = varobs;

% Test the number of variables in the database.
if dataset_.info.nvobs-size(rawdata,2)
    disp(' ')
    disp(['Declared number of observed variables = ' int2str(dataset.info.nvobs)])
    disp(['Number of variables in the database   = ' int2str(size(rawdata,2))])
    disp(' ')
    error(['Estimation can''t take place because the declared number of observed' ...
           'variables doesn''t match the number of variables in the database.'])
end

if size(rawdata,1)~=dataset_.info.ntobs
   fprintf('Restricting the sample to observations %d to %d. Using in total %d observations. \n',first,dataset_.info.ntobs,dataset_.info.ntobs-first+1)
end
rawdata = rawdata(first:(first+dataset_.info.ntobs-1),:);

% Take the log (or anything else) of the variables if needed
if isempty(transformation)
    dataset_.rawdata = rawdata;
else
    dataset_.rawdata = arrayfun(transformation,rawdata);
end

% Test if the observations are real numbers.
if ~isreal(dataset_.rawdata)
    error('Estimation:: There are complex values in the data! Probably  a wrong (log) transformation...')
end

% Test for missing observations.
dataset_.missing.state = any(any(isnan(dataset_.rawdata)));
if dataset_.missing.state
    [i,n,s,j] = describe_missing_data(dataset_.rawdata);
    dataset_.missing.aindex = i;
    dataset_.missing.vindex = j;
    dataset_.missing.number_of_observations = n;
    dataset_.missing.no_more_missing_observations = s;
else
    dataset_.missing.aindex = num2cell(repmat(1:dataset_.info.nvobs,dataset_.info.ntobs,1)',1);
    dataset_.missing.vindex = [];
    dataset_.missing.number_of_observations = [];
    dataset_.missing.no_more_missing_observations = 1;
end

% Compute the empirical mean of the observed variables..
dataset_.descriptive.mean = nanmean(dataset_.rawdata);

% Prefilter the data if needed.
if prefilter == 1
    dataset_.data = transpose(bsxfun(@minus,dataset_.rawdata,dataset_.descriptive.mean));
else
    dataset_.data = transpose(dataset_.rawdata);
end
