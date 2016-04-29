function DynareOutput = simul_backward_nonlinear_model(initial_conditions, sample_size, DynareOptions, DynareModel, DynareOutput, innovations)

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
    DynareOutput.exo_simul = zeros(sample_size,number_of_shocks);
    DynareOutput.exo_simul(:,positive_var_indx) = DynareOutput.bnlms.shocks;
    DynareOutput.exo_simul = [zeros(1,number_of_shocks); DynareOutput.exo_simul];
else
    DynareOutput.exo_simul = innovations;
end

% Get usefull vector of indices.
ny0 = nnz(DynareModel.lead_lag_incidence(2,:));
ny1 = nnz(DynareModel.lead_lag_incidence(1,:));
iy1 = find(DynareModel.lead_lag_incidence(1,:)>0);
idx = 1:DynareModel.endo_nbr;
jdx = idx+ny1;
hdx = 1:ny1;

% Get the name of the dynamic model routine.
model_dynamic = str2func([DynareModel.fname,'_dynamic']);

% initialization of vector y.
y = NaN(length(idx)+ny1,1);

% initialization of the returned simulations.
DynareOutput.endo_simul = NaN(DynareModel.endo_nbr,sample_size+1);
if isempty(initial_conditions)
    DynareOutput.endo_simul(:,1) = DynareOutput.steady_state;
else
    DynareOutput.endo_simul(:,1) = initial_conditions;
end
Y = DynareOutput.endo_simul;

% Simulations (call a Newton-like algorithm for each period).
for it = 2:sample_size+1
    y(jdx) = Y(:,it-1);                       % A good guess for the initial conditions is the previous values for the endogenous variables.
    y(hdx) = y(jdx(iy1));                     % Set lagged variables.
    z = trust_region(model_dynamic, y, idx, jdx, 1, DynareOptions.gstep, ...
                    DynareOptions.solve_tolf,DynareOptions.solve_tolx, ...
                    DynareOptions.simul.maxit,DynareOptions.debug, ...
                    DynareOutput.exo_simul, DynareModel.params, ...
                    DynareOutput.steady_state, it);
    Y(:,it) = z(jdx);
end

DynareOutput.endo_simul = Y;