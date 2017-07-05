function [nvar,vartan,NumberOfConditionalDecompFiles] = ...
    dsge_simulated_theoretical_conditional_variance_decomposition(SampleSize,Steps,M_,options_,oo_,type)
% function [nvar,vartan,NumberOfConditionalDecompFiles] = ...
%     dsge_simulated_theoretical_conditional_variance_decomposition(SampleSize,Steps,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the conditional variance
% decomposition of the endogenous variables (or a subset of the endogenous variables).
%
% INPUTS
%   SampleSize   [integer]       scalar, number of simulations.
%   Steps        [integers]      horizons at which to conduct decomposition
%   M_           [structure]     Dynare structure describing the model.
%   options_     [structure]     Dynare structure defining global options.
%   oo_          [structure]     Dynare structure where the results are saved.
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   nvar                             [integer]  nvar is the number of stationary variables.
%   vartan                           [char]     array of characters (with nvar rows).
%   NumberOfConditionalDecompFiles   [integer]  scalar, number of prior or posterior data files (for covariance).

% Copyright (C) 2009-2015 Dynare Team
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


% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
    posterior = 1;
elseif strcmpi(type,'prior')
    DrawsFiles = dir([M_.dname '/prior/draws/' type '_draws*' ]);
    CheckPath('prior/moments',M_.dname);
    posterior = 0;
else
    disp('dsge_simulated_theoretical_conditional_variance_decomposition:: Unknown type!')
    error()
end

%delete old stale files before creating new ones
if posterior
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_PosteriorConditionalVarianceDecomposition*'])
else
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_PriorConditionalVarianceDecomposition*'])
end

% Set varlist (vartan)
if ~posterior
    if isfield(options_,'varlist')
        temp = options_.varlist;
    end
    options_.varlist = options_.prior_analysis_endo_var_list;
end
[ivar,vartan ] = get_variables_list(options_, M_);
if ~posterior
    if exist('temp','var')
        options_.varlist = temp;
    end
end
nvar = length(ivar);

% Set the size of the auto-correlation function to zero.
nar = options_.ar;
options_.ar = 0;

NumberOfDrawsFiles = rows(DrawsFiles);
NumberOfSavedElementsPerSimulation = nvar*M_.exo_nbr*length(Steps);
MaXNumberOfConditionalDecompLines = ceil(options_.MaxNumberOfBytes/NumberOfSavedElementsPerSimulation/8);

if SampleSize<=MaXNumberOfConditionalDecompLines
    Conditional_decomposition_array = zeros(nvar,length(Steps),M_.exo_nbr,SampleSize);
    NumberOfConditionalDecompFiles = 1;
else
    Conditional_decomposition_array = zeros(nvar,length(Steps),M_.exo_nbr,MaXNumberOfConditionalDecompLines);
    NumberOfLinesInTheLastConditionalDecompFile = mod(SampleSize,MaXNumberOfConditionalDecompLines);
    NumberOfConditionalDecompFiles = ceil(SampleSize/MaXNumberOfConditionalDecompLines);
end

NumberOfConditionalDecompLines = size(Conditional_decomposition_array,4);
ConditionalDecompFileNumber = 0;

StateSpaceModel.number_of_state_equations = M_.endo_nbr;
StateSpaceModel.number_of_state_innovations = M_.exo_nbr;

first_call = 1;

linea = 0;
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
            M_=set_parameters_locally(M_,pdraws{linee,1});% Needed to update the covariance matrix of the state innovations.
            dr = pdraws{linee,2};
        else
            M_=set_parameters_locally(M_,pdraws{linee,1});
            [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        end
        if first_call
            endo_nbr = M_.endo_nbr;
            nstatic = M_.nstatic;
            nspred = M_.nspred;
            iv = (1:endo_nbr)';
            ic = [ nstatic+(1:nspred) endo_nbr+(1:size(dr.ghx,2)-nspred) ]';
            StateSpaceModel.number_of_state_equations = M_.endo_nbr;
            StateSpaceModel.number_of_state_innovations = M_.exo_nbr;
            StateSpaceModel.sigma_e_is_diagonal = M_.sigma_e_is_diagonal;
            StateSpaceModel.order_var = dr.order_var;
            first_call = 0;
            clear('endo_nbr','nstatic','nspred','k');
        end
        [StateSpaceModel.transition_matrix,StateSpaceModel.impulse_matrix] = kalman_transition_matrix(dr,iv,ic,M_.exo_nbr);
        StateSpaceModel.state_innovations_covariance_matrix = M_.Sigma_e;
        clear('dr');
        Conditional_decomposition_array(:,:,:,linea) = conditional_variance_decomposition(StateSpaceModel, Steps, ivar);
        if linea == NumberOfConditionalDecompLines
            ConditionalDecompFileNumber = ConditionalDecompFileNumber + 1;
            linea = 0;
            if posterior
                save([M_.dname '/metropolis/' M_.fname '_PosteriorConditionalVarianceDecomposition' int2str(ConditionalDecompFileNumber) '.mat' ], ...
                     'Conditional_decomposition_array');
            else
                save([M_.dname '/prior/moments/' M_.fname '_PriorConditionalVarianceDecomposition' int2str(ConditionalDecompFileNumber) '.mat' ], ...
                     'Conditional_decomposition_array');
            end
            if (ConditionalDecompFileNumber==NumberOfConditionalDecompFiles-1)% Prepare last round.
                Conditional_decomposition_array = zeros(nvar, length(Steps),M_.exo_nbr,NumberOfLinesInTheLastConditionalDecompFile) ;
                NumberOfConditionalDecompLines = NumberOfLinesInTheLastConditionalDecompFile;
            elseif ConditionalDecompFileNumber<NumberOfConditionalDecompFiles-1
                Conditional_decomposition_array = zeros(nvar,length(Steps),M_.exo_nbr,MaXNumberOfConditionalDecompLines);
            else
                clear('Conditional_decomposition_array');
            end
        end
    end
end

options_.ar = nar;
