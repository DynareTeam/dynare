function [nvar,vartan,CorrFileNumber] = dsge_simulated_theoretical_correlation(SampleSize,nar,M_,options_,oo_,type)
% function [nvar,vartan,CorrFileNumber] = dsge_simulated_theoretical_correlation(SampleSize,nar,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the endogenous
% variables' second order moments.
%
% INPUTS
%   SampleSize   [integer]          scalar, number of simulations.
%   nar          [integer]          maximum number of autocorrelations to
%                                   consider
%   M_           [structure]        Dynare structure describing the model.
%   options_     [structure]        Dynare structure defining global options
%   oo_          [structure]        Dynare structure where the results are saved.
%   type         [string]           'prior' or 'posterior'
%
% OUTPUTS
%   nvar           [integer]        nvar is the number of stationary variables.
%   vartan         [char]           array of characters (with nvar rows).
%   CorrFileNumber [integer]        scalar, number of prior or posterior data files (for correlation).

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
    disp('dsge_simulated_theoretical_correlation:: Unknown type!');
    error()
end
NumberOfDrawsFiles = length(DrawsFiles);

%delete old stale files before creating new ones
if posterior
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations*']);
else
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_PriorCorrelations*']);
end

% Set varlist (vartan)
if ~posterior
    if isfield(options_,'varlist')
        temp = options_.varlist;
    end
    options_.varlist = options_.prior_analysis_endo_var_list;
end
[ivar,vartan, options_] = get_variables_list(options_, M_);
if ~posterior
    if exist('temp','var')
        options_.varlist = temp;
    end
end
nvar = length(ivar);

% Set the size of the auto-correlation function to nar.
oldnar = options_.ar;
options_.ar = nar;

% Number of lines in posterior data files.
MaXNumberOfCorrLines = ceil(options_.MaxNumberOfBytes/(nvar*nvar*nar)/8);

if SampleSize<=MaXNumberOfCorrLines
    Correlation_array = zeros(SampleSize,nvar,nvar,nar);
    NumberOfCorrFiles = 1;
else
    Correlation_array = zeros(MaXNumberOfCorrLines,nvar,nvar,nar);
    NumberOfLinesInTheLastCorrFile = mod(SampleSize,MaXNumberOfCorrLines);
    NumberOfCorrFiles = ceil(SampleSize/MaXNumberOfCorrLines);
end

NumberOfCorrLines = rows(Correlation_array);
CorrFileNumber = 1;

% Compute 2nd order moments and save them in *_[Posterior, Prior]Correlations* files
linea = 0;
for file = 1:NumberOfDrawsFiles
    if posterior
        load([M_.dname '/metropolis/' DrawsFiles(file).name ]);
    else
        load([M_.dname '/prior/draws/' DrawsFiles(file).name]);
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
        for i=1:nar
            Correlation_array(linea,:,:,i) = tmp{i+1};
        end
        if linea == NumberOfCorrLines
            if posterior
                save([ M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations' int2str(CorrFileNumber) '.mat' ],'Correlation_array');
            else
                save([ M_.dname '/prior/moments/' M_.fname '_PriorCorrelations' int2str(CorrFileNumber) '.mat' ],'Correlation_array');
            end
            CorrFileNumber = CorrFileNumber + 1;
            linea = 0;
            test = CorrFileNumber-NumberOfCorrFiles;
            if ~test% Prepare the last round...
                Correlation_array = zeros(NumberOfLinesInTheLastCorrFile,nvar,nvar,nar);
                NumberOfCorrLines = NumberOfLinesInTheLastCorrFile;
                CorrFileNumber = CorrFileNumber - 1;
            elseif test<0
                Correlation_array = zeros(MaXNumberOfCorrLines,nvar,nvar,nar);
            else
                clear('Correlation_array');
            end
        end
    end
end

options_.ar = oldnar;