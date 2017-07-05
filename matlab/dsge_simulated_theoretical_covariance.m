function [nvar,vartan,CovarFileNumber] = dsge_simulated_theoretical_covariance(SampleSize,M_,options_,oo_,type)
% function [nvar,vartan,CovarFileNumber] = dsge_simulated_theoretical_covariance(SampleSize,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the endogenous
% variables second order moments.
%
% INPUTS
%   SampleSize   [integer]       scalar, number of simulations.
%   M_           [structure]     Dynare structure describing the model.
%   options_     [structure]     Dynare structure defining global options.
%   oo_          [structure]     Dynare structure where the results are saved.
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   nvar              [integer]  nvar is the number of stationary variables.
%   vartan            [char]     array of characters (with nvar rows).
%   CovarFileNumber   [integer]  scalar, number of prior or posterior data files (for covariance).

% Copyright (C) 2007-2017 Dynare Team
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

nodecomposition = 1;

% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
    posterior = 1;
elseif strcmpi(type,'prior')
    DrawsFiles = dir([M_.dname '/prior/draws/' type '_draws*' ]);
    CheckPath('prior/moments',M_.dname);
    posterior = 0;
else
    disp('dsge_simulated_theoretical_covariance:: Unknown type!')
    error();
end
NumberOfDrawsFiles = length(DrawsFiles);

%delete old stale files before creating new ones
if posterior
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments*'])
else
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_Prior2ndOrderMoments*'])
end

% Set varlist (vartan)
if ~posterior
    if isfield(options_,'varlist')
        temp = options_.varlist;
    end
    options_.varlist = options_.prior_analysis_endo_var_list;
end
[ivar,vartan] = get_variables_list(options_,M_);
if ~posterior
    if exist('temp','var')
        options_.varlist = temp;
    end
end
nvar = length(ivar);

% Set the size of the auto-correlation function to zero.
nar = options_.ar;
options_.ar = 0;

% Number of lines in posterior data files.
MaXNumberOfCovarLines = ceil(options_.MaxNumberOfBytes/(nvar*(nvar+1)/2)/8);

if SampleSize<=MaXNumberOfCovarLines
    Covariance_matrix = zeros(SampleSize,nvar*(nvar+1)/2);
    NumberOfCovarFiles = 1;
else
    Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
    NumberOfLinesInTheLastCovarFile = mod(SampleSize,MaXNumberOfCovarLines);
    NumberOfCovarFiles = ceil(SampleSize/MaXNumberOfCovarLines);
end

NumberOfCovarLines = rows(Covariance_matrix);
CovarFileNumber = 1;

% Compute 2nd order moments and save them in *_[Posterior, Prior]2ndOrderMoments* files
linea = 0;
for file = 1:NumberOfDrawsFiles
    if posterior
        load([M_.dname '/metropolis/' DrawsFiles(file).name ],'pdraws');
    else
        load([M_.dname '/prior/draws/' DrawsFiles(file).name ],'pdraws');
    end
    NumberOfDraws = rows(pdraws);
    isdrsaved = columns(pdraws)-1;
    for linee = 1:NumberOfDraws
        linea = linea+1;
        if isdrsaved
            dr = pdraws{linee,2};
        else
            M_=set_parameters_locally(M_,pdraws{linee,1});
            [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        end
        tmp = th_autocovariances(dr,ivar,M_,options_,nodecomposition);
        for i=1:nvar
            for j=i:nvar
                Covariance_matrix(linea,symmetric_matrix_index(i,j,nvar)) = tmp{1}(i,j);
            end
        end
        if linea == NumberOfCovarLines
            if posterior
                save([ M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments' int2str(CovarFileNumber) '.mat' ],'Covariance_matrix');
            else
                save([ M_.dname '/prior/moments/' M_.fname '_Prior2ndOrderMoments' int2str(CovarFileNumber) '.mat' ],'Covariance_matrix');
            end
            CovarFileNumber = CovarFileNumber + 1;
            linea = 0;
            test = CovarFileNumber-NumberOfCovarFiles;
            if ~test% Prepare the last round...
                Covariance_matrix = zeros(NumberOfLinesInTheLastCovarFile,nvar*(nvar+1)/2);
                NumberOfCovarLines = NumberOfLinesInTheLastCovarFile;
            elseif test<0
                Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
            else
                clear('Covariance_matrix');
            end
        end
    end
end

options_.ar = nar;