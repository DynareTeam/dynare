function [nvar,vartan,NumberOfDecompFiles] = ...
    dsge_simulated_theoretical_variance_decomposition(SampleSize,M_,options_,oo_,type)
% function [nvar,vartan,NumberOfDecompFiles] = ...
%     dsge_simulated_theoretical_variance_decomposition(SampleSize,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the variance
% decomposition of the observed endogenous variables.
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

nodecomposition = 0;

% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
    posterior = 1;
elseif strcmpi(type,'prior')
    DrawsFiles = dir([M_.dname '/prior/draws/' type '_draws*' ]);
    CheckPath('prior/moments',M_.dname);
    posterior = 0;
else
    disp('dsge_simulated_theoretical_variance_decomposition:: Unknown type!')
    error()
end

%delete old stale files before creating new ones
if posterior
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_PosteriorVarianceDecomposition*']);
else
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_PosteriorVarianceDecomposition*']);
end

% Set varlist (vartan)
if ~posterior
    if isfield(options_,'varlist')
        temp = options_.varlist;
    end
    options_.varlist = options_.prior_analysis_endo_var_list;
end
[ivar,vartan,options_] = get_variables_list(options_,M_);
if ~posterior
    if exist('temp','var')
        options_.varlist = temp;
    end
end
nvar = length(ivar);

% Set the size of the auto-correlation function to zero.
nar = options_.ar;
options_.ar = 0;



nexo = M_.exo_nbr;

NumberOfDrawsFiles = rows(DrawsFiles);
NumberOfSavedElementsPerSimulation = nvar*(nexo+1);
MaXNumberOfDecompLines = ceil(options_.MaxNumberOfBytes/NumberOfSavedElementsPerSimulation/8);

if SampleSize<=MaXNumberOfDecompLines
    Decomposition_array = zeros(SampleSize,nvar*nexo);
    NumberOfDecompFiles = 1;
else
    Decomposition_array = zeros(MaXNumberOfDecompLines,nvar*nexo);
    NumberOfLinesInTheLastDecompFile = mod(SampleSize,MaXNumberOfDecompLines);
    NumberOfDecompFiles = ceil(SampleSize/MaXNumberOfDecompLines);
end

NumberOfDecompLines = rows(Decomposition_array);
DecompFileNumber = 1;

% Compute total variances (covariances are not saved) and variances
% implied by each structural shock.
linea = 0;
only_non_stationary_vars=0;
for file = 1:NumberOfDrawsFiles
    if posterior
        load([M_.dname '/metropolis/' DrawsFiles(file).name ]);
    else
        load([M_.dname '/prior/draws/' DrawsFiles(file).name ]);
    end
    isdrsaved = columns(pdraws)-1;
    NumberOfDraws = rows(pdraws);
    for linee = 1:NumberOfDraws
        linea = linea+1;
        if isdrsaved
            dr = pdraws{linee,2};
        else
            M_=set_parameters_locally(M_,pdraws{linee,1});
            [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        end
        if file==1 && linee==1
            [tmp, stationary_vars] = th_autocovariances(dr,ivar,M_,options_,nodecomposition);
            if isempty(stationary_vars)
                fprintf('\ndsge_simulated_theoretical_variance_decomposition:: All requested endogenous variables have a unit root and thus infinite variance.\n')
                fprintf('dsge_simulated_theoretical_variance_decomposition:: No decomposition is performed.\n')
                only_non_stationary_vars=1;
            end
        end
        if only_non_stationary_vars
            for i=1:nvar
                for j=1:nexo
                    Decomposition_array(linea,(i-1)*nexo+j) = NaN;
                end
            end
        else
            tmp = th_autocovariances(dr,ivar,M_,options_,nodecomposition);
            for i=1:nvar
                for j=1:nexo
                    Decomposition_array(linea,(i-1)*nexo+j) = tmp{2}(i,j);
                end
            end
        end
        if linea == NumberOfDecompLines
            if posterior
                save([M_.dname '/metropolis/' M_.fname '_PosteriorVarianceDecomposition' int2str(DecompFileNumber) '.mat' ],'Decomposition_array');
            else
                save([M_.dname '/prior/moments/' M_.fname '_PriorVarianceDecomposition' int2str(DecompFileNumber) '.mat' ],'Decomposition_array');
            end
            DecompFileNumber = DecompFileNumber + 1;
            linea = 0;
            test = DecompFileNumber-NumberOfDecompFiles;
            if ~test% Prepare the last round...
                Decomposition_array = zeros(NumberOfLinesInTheLastDecompFile,nvar*nexo);
                NumberOfDecompLines = NumberOfLinesInTheLastDecompFile;
            elseif test<0
                Decomposition_array = zeros(MaXNumberOfDecompLines,nvar*nexo);
            else
                clear('Decomposition_array');
            end
        end
    end
end

options_.ar = nar;