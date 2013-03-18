function record=independent_metropolis_hastings(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin)

% Independent Metropolis-Hastings algorithm. 
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o varargin              list of argument following mh_bounds
%  
% OUTPUTS 
%   o record     [struct]   structure describing the iterations
%
% ALGORITHM 
%   Metropolis-Hastings.       
%
% SPECIAL REQUIREMENTS
%   None.
%
% PARALLEL CONTEXT
% See the comment in random_walk_metropolis_hastings.m funtion.

% Copyright (C) 2006-2011 Dynare Team
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

global M_ options_ bayestopt_ estim_params_ oo_
%%%%
%%%% Initialization of the independent metropolis-hastings chains.
%%%%
[ ix2, ilogpo2, ModelName, MhDirectoryName, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = ...
    metropolis_hastings_initialization(TargetFun, xparam1, vv, mh_bounds, varargin{:});

xparam1 = transpose(xparam1);      
InitSizeArray = min([repmat(MAX_nruns,nblck,1) fline+nruns-1],[],2);

load([MhDirectoryName '/' ModelName '_mh_history.mat'],'record');

%The mandatory variables for local/remote parallel computing are stored in localVars struct.

localVars =   struct('TargetFun', TargetFun, ...
                     'ProposalFun', ProposalFun, ...
                     'xparam1', xparam1, ...
                     'vv', vv, ...
                     'mh_bounds', mh_bounds, ...
                     'ix2', ix2, ...
                     'ilogpo2', ilogpo2, ...
                     'ModelName', ModelName, ...
                     'fline', fline, ...
                     'npar', npar, ...
                     'nruns', nruns, ...
                     'NewFile', NewFile, ...
                     'MAX_nruns', MAX_nruns, ...
                     'd', d);
localVars.InitSizeArray=InitSizeArray;
localVars.record=record;
localVars.varargin=varargin;

% Like a sequential execution!
if isnumeric(options_.parallel),
    fout = independent_metropolis_hastings_core(localVars, fblck, nblck, 0);
    record = fout.record;
    % Parallel execution.    
else
    % global variables for parallel routines
    globalVars = struct('M_',M_, ...
                        'options_', options_, ...
                        'bayestopt_', bayestopt_, ...
                        'estim_params_', estim_params_, ...
                        'oo_', oo_);
    
    % which files have to be copied to run remotely
    NamFileInput(1,:) = {'',[ModelName '_static.m']};
    NamFileInput(2,:) = {'',[ModelName '_dynamic.m']};
    if options_.steadystate_flag,
        NamFileInput(length(NamFileInput)+1,:)={'',[ModelName '_steadystate.m']};
    end
    if (options_.load_mh_file~=0) && any(fline>1) ,
        NamFileInput(length(NamFileInput)+1,:)={[M_.dname '/metropolis/'],[ModelName '_mh' int2str(NewFile(1)) '_blck*.mat']};
    end
    if exist([ModelName '_optimal_mh_scale_parameter.mat'])
        NamFileInput(length(NamFileInput)+1,:)={'',[ModelName '_optimal_mh_scale_parameter.mat']};
    end
    
    % from where to get back results
    %     NamFileOutput(1,:) = {[M_.dname,'/metropolis/'],'*.*'};
    
    [fout, nBlockPerCPU, totCPU] = masterParallel(options_.parallel, fblck, nblck,NamFileInput,'independent_metropolis_hastings_core', localVars, globalVars, options_.parallel_info);
    for j=1:totCPU,
        offset = sum(nBlockPerCPU(1:j-1))+fblck-1;
        record.LastLogPost(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.LastLogPost(offset+1:sum(nBlockPerCPU(1:j)));
        record.LastParameters(offset+1:sum(nBlockPerCPU(1:j)),:)=fout(j).record.LastParameters(offset+1:sum(nBlockPerCPU(1:j)),:);
        record.AcceptationRates(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.AcceptationRates(offset+1:sum(nBlockPerCPU(1:j)));
        record.Seeds(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.Seeds(offset+1:sum(nBlockPerCPU(1:j)));
    end

end

irun = fout(1).irun;
NewFile = fout(1).NewFile;

% record.Seeds.Normal = randn('state');
% record.Seeds.Unifor = rand('state');
save([MhDirectoryName '/' ModelName '_mh_history.mat'],'record');
disp(['MH: Number of mh files                   : ' int2str(NewFile(1)) ' per block.'])
disp(['MH: Total number of generated files      : ' int2str(NewFile(1)*nblck) '.'])
disp(['MH: Total number of iterations           : ' int2str((NewFile(1)-1)*MAX_nruns+irun-1) '.'])
disp('MH: average acceptation rate per chain   : ')
disp(record.AcceptationRates);
disp(' ')