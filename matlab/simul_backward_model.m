function DynareOutput = simul_backward_linear_model(initial_conditions, sample_size, DynareOptions, DynareModel, DynareOutput, innovations)

%@info:
%! @deftypefn {Function File} {@var{DynareOutput} =} simul_backward_nonlinear_model (@var{sample_size},@var{DynareOptions}, @var{DynareModel}, @var{DynareOutput})
%! @anchor{@simul_backward_nonlinear_model}
%! @sp 1
%! Simulates a stochastic non linear backward looking model with arbitrary precision (a deterministic solver is used).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item sample_size
%! Scalar integer, size of the sample to be generated.
%! @item DynareOptions
%! Matlab/Octave structure (Options used by Dynare).
%! @item DynareDynareModel
%! Matlab/Octave structure (Description of the model).
%! @item DynareOutput
%! Matlab/Octave structure (Results reported by Dynare).
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item DynareOutput
%! Matlab/Octave structure (Results reported by Dynare).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{dynTime}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2012-2016 Dynare Team
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

if DynareModel.maximum_lead
    error(['simul_backward_nonlinear_model:: The specified model is not backward looking!'])
end

if nargin<6
    % Set the covariance matrix of the structural innovations.
    variances = diag(DynareModel.Sigma_e);
    number_of_shocks = length(DynareModel.Sigma_e);
    positive_var_indx = find(variances>0);
    effective_number_of_shocks = length(positive_var_indx);
    covariance_matrix = DynareModel.Sigma_e(positive_var_indx,positive_var_indx);
    covariance_matrix_upper_cholesky = chol(covariance_matrix);

    % Set seed to its default state.
    if DynareOptions.bnlms.set_dynare_seed_to_default
        set_dynare_seed('default');
    end

    % Simulate structural innovations.
    switch DynareOptions.bnlms.innovation_distribution
      case 'gaussian'
        DynareOutput.bnlms.shocks = randn(sample_size,effective_number_of_shocks)*covariance_matrix_upper_cholesky;
      otherwise
        error(['simul_backward_nonlinear_model:: ' DynareOption.bnlms.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end

    % Put the simulated innovations in DynareOutput.exo_simul.
    DynareOutput.exo_simul = zeros(sample_size+1,number_of_shocks);
    DynareOutput.exo_simul(2:end,positive_var_indx) = DynareOutput.bnlms.shocks;
else
    number_of_shocks = size(innovations,2);
    DynareOutput.exo_simul = innovations;
end

if DynareOptions.linear
    DynareOutput = simul_backward_linear_model(initial_conditions, sample_size, DynareOptions, ...
                                DynareModel, DynareOutput, innovations);
else
    DynareOutput = simul_backward_nonlinear_model(initial_conditions, sample_size, DynareOptions, ...
                                DynareModel, DynareOutput, innovations);
end    
    