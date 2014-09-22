function [DynareDataset, DatasetInfo, newdatainterface] = makedataset(DynareOptions, initialconditions, gsa_flag)

% Initialize a dataset as a dseries object.
%
%
% INPUTS
% ======
%
%     DynareOptions [struct] Structure of options built by Dynare's preprocessor.
%
%
% OUTPUTS
% =======
%
%     DynareDataset [dseries]  The dataset.
%     DatasetInfo   [struct]   Various informations about the dataset (descriptive statistics and missing observations).
%
% EXAMPLE
% =======
%
%     [dataset_, dataset_info] = makedataset(options_) ;
%
%
% See also dynare_estimation_init

if nargin<3
    gsa_flag = 0;
end

if nargin<2 || isempty(initialconditions)
    % If a the sample is to be used for the estimation of a VAR or DSGE-VAR model
    % the second argument must be a strictly positive integer (the number of lags).
    initialconditions = 0;
end

if isempty(DynareOptions.datafile) && isempty(DynareOptions.dataset.file) && isempty(DynareOptions.dataset.series)
    if gsa_flag
        DynareDataset = dseries();
        DatasetInfo = [];
        newdatainterface=0;
        return
    else
        error('makedataset: datafile option is missing!')
    end
end

if isempty(DynareOptions.datafile) && ~isempty(DynareOptions.dataset.file)
    datafile = DynareOptions.dataset.file;
    newdatainterface = 1;
elseif isempty(DynareOptions.datafile) && ~isempty(DynareOptions.dataset.series)
    try
        dseriesobjectforuserdataset = evalin('base', DynareOptions.dataset.series);
    catch
        error(sprintf('makedataset: %s is unknown!', DynareOptions.dataset.series))
    end
    if ~isdseries(dseriesobjectforuserdataset)
        error(sprintf('makedataset: %s has to be a dseries object!', DynareOptions.dataset.series))
    end
    datafile = [];
    newdatainterface = 1;
elseif ~isempty(DynareOptions.datafile) && isempty(DynareOptions.dataset.file)
    datafile = DynareOptions.datafile;
    newdatainterface = 0;
elseif isempty(DynareOptions.datafile) && ~isempty(DynareOptions.dataset.file)
    error('makedataset: You cannot simultaneously use the data command and the datafile option (in the estimation command)!')
else
    error('makedataset: You have to specify the datafile!')
end

% Check extension.
if ~isempty(datafile)
    allowed_extensions = {'m','mat','csv','xls','xlsx'};
    datafile_extension = get_file_extension(datafile);
    if isempty(datafile_extension)
        available_extensions = {}; j = 1;
        for i=1:length(allowed_extensions)
            if exist([datafile '.' allowed_extensions{i}])
                available_extensions(j) = {allowed_extensions{i}};
                j = j+1;
            end
        end
        if isempty(available_extensions)
            error(['makedataset: I can''t find a datafile (with allowed extension m, mat, csv, xls or xlsx)!'])
        end
        if length(available_extensions)>1
            error(sprintf(['makedataset: You did not specify an extension for the datafile, but more than one candidate ' ...
                           'is available in the designated folder!\nPlease, add an extension to the datafile ' ...
                           '(m, mat, csv, xls or xlsx are permitted extensions).']));
        end
        datafile = [datafile '.' available_extensions{1}];
    end
end

% Load the data in a dseries object.
if ~isempty(datafile)
    DynareDataset = dseries(datafile);
else
    DynareDataset = dseriesobjectforuserdataset;
    clear('dseriesobjectforuserdataset');
end

% Select a subset of the variables.
DynareDataset = DynareDataset{DynareOptions.varobs{:}};

% Apply log function if needed.
if DynareOptions.loglinear && ~DynareOptions.logdata
    DynareDataset = DynareDataset.log();
end

% Test if an initial period (different from its default value) is explicitely defined in the datafile.
if isequal(DynareDataset.init, dates(1,1))
    dataset_default_initial_period = 1;
else
    dataset_default_initial_period = 0;
end

%  Test if an initial period (different from its default value) is explicitely defined in the mod file with the set_time command.
if ~isdates(DynareOptions.initial_period) && isnan(DynareOptions.initial_period)
    set_time_default_initial_period = 1;
else
    set_time_default_initial_period = 0;
end

if ~set_time_default_initial_period && dataset_default_initial_period
    % Overwrite the initial period in dataset (it was set to default).
    % Note that the updates of freq and time members are auto-magically
    % done by dseries::subsasgn overloaded method.
    DynareDataset.init = DynareOptions.initial_period;
end

if set_time_default_initial_period && ~dataset_default_initial_period
    % Overwrite the global initial period defined by set_time (it was set to default).
    DynareOptions.initial_period = DynareDataset.init;
end

if ~set_time_default_initial_period && ~dataset_default_initial_period
    % Check if dataset.init and options_.initial_period are identical.
    if DynareOptions.initial_period<DynareDataset.init
        error('makedataset: The date as defined by the set_time command is not consistent with the initial period in the database!')
    end
end

% Set firstobs, lastobs and nobs
if newdatainterface
    if isempty(DynareOptions.dataset.firstobs)
        % first_obs option was not used in the data command.
        firstobs = DynareDataset.init;
    else
        firstobs = DynareOptions.dataset.firstobs;
    end
    if isnan(DynareOptions.dataset.nobs)
        % nobs option was not used in the data command.
        if isempty(DynareOptions.dataset.lastobs)
            % last_obs option was not used in the data command.
            nobs = DynareDataset.nobs;
            lastobs = DynareDataset.dates(end);
        else
            lastobs = DynareOptions.dataset.lastobs;
            nobs = lastobs-firstobs+1;
        end
    else
        nobs = DynareOptions.dataset.nobs;
        if isempty(DynareOptions.dataset.lastobs)
            % last_obs option was not used in the data command.
            lastobs = firstobs+(nobs-1);
        else
            % last_obs and nobs were used in the data command. Check that they are consistent (with firstobs).
            if ~isequal(lastobs,firstobs+(nobs-1))
                error(sprintf('makedataset: Options last_obs (%s), first_obs (%s) and nobs (%s) are not consistent!',char(lastobs),char(firstobs),num2str(nobs)));
            end
        end
    end
else
    if isnan(DynareOptions.first_obs)
        firstobs = DynareDataset.init;
    else
        firstobs = DynareDataset.dates(DynareOptions.first_obs);
    end
    if isnan(DynareOptions.nobs)
        lastobs = DynareDataset.dates(end);
        nobs = lastobs-firstobs+1;
    else
        nobs = DynareOptions.nobs;
        lastobs = firstobs+(nobs-1);
    end
end

% Add initial conditions if needed
FIRSTOBS = firstobs-initialconditions;

% Check that firstobs belongs to DynareDataset.dates
if firstobs<DynareDataset.init
    error(sprintf('makedataset: first_obs (%s) cannot be less than the first date in the dataset (%s)!',char(firstobs),char(DynareDataset.init)))
end

% Check that FIRSTOBS belongs to DynareDataset.dates
if initialconditions && FIRSTOBS<DynareDataset.init
    error(sprintf('makedataset: first_obs (%s) - %i cannot be less than the first date in the dataset (%s)!\nReduce the number of lags in the VAR model or increase the value of first_obs\nto at least first_obs=%i.', char(firstobs), initialconditions, char(DynareDataset.init),initialconditions+1));
end

% Check that lastobs belongs to DynareDataset.dates...
if newdatainterface
    if lastobs>DynareDataset.dates(end)
        error(sprintf('makedataset: last_obs (%s) cannot be greater than the last date in the dataset (%s)!',char(lastobs),char(DynareDataset.dates(end))))
    end
else
    % ...  or check that nobs is smaller than the number of observations in DynareDataset.
    if nobs>DynareDataset.nobs
        error(sprintf('makedataset: nobs (%s) cannot be greater than the last date in the dataset (%s)!', num2str(nobs), num2str(DynareDataset.nobs)))
    end
end

% Select a subsample.
DynareDataset = DynareDataset(FIRSTOBS:lastobs);

% Initialize DatasetInfo structure.
DatasetInfo = struct('missing', struct('state', NaN, 'aindex', [], 'vindex', [], 'number_of_observations', NaN, 'no_more_missing_observations', NaN), ...
                     'descriptive', struct('mean', [], 'covariance', [], 'correlation', [], 'autocovariance', []));

% Fill DatasetInfo.missing if some observations are missing
DatasetInfo.missing.state = isanynan(DynareDataset.data);
if DatasetInfo.missing.state
    [DatasetInfo.missing.aindex, DatasetInfo.missing.number_of_observations, DatasetInfo.missing.no_more_missing_observations, DatasetInfo.missing.vindex] = ...
        describe_missing_data(DynareDataset.data);
else
    DatasetInfo.missing.aindex = num2cell(transpose(repmat(1:DynareDataset.vobs,DynareDataset.nobs,1)),1);
    DatasetInfo.missing.no_more_missing_observations = 1;
end

% Compute the empirical mean of the observed variables.
DatasetInfo.descriptive.mean = nanmean(DynareDataset.data);

% Compute the empirical covariance matrix of the observed variables.
DatasetInfo.descriptive.covariance = nancovariance(DynareDataset.data);

% Compute the empirical correlation matrix of the observed variables.
normalization_matrix = diag(1./sqrt(diag(DatasetInfo.descriptive.covariance)));
DatasetInfo.descriptive.correlation = normalization_matrix*DatasetInfo.descriptive.covariance*normalization_matrix;

% Compute autocorrelation function.
DatasetInfo.descriptive.autocovariance = nanautocovariance(DynareDataset.data, DynareOptions.ar);

% Save raw data.
DatasetInfo.rawdata = DynareDataset.data;

% Prefilter the data if needed (remove the mean).
if isequal(DynareOptions.prefilter, 1)
    DynareDataset = DynareDataset.detrend();
end
