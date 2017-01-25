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

number_of_shocks = size(DynareOutput.exo_simul,2);
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

% get coefficients
[cst,jacob] = model_dynamic(zeros(DynareModel.endo_nbr+ny1,1), ...
                            zeros(2,size(DynareOutput.exo_simul, 2)), ...
                            DynareModel.params, ...
                            DynareOutput.steadystate,2);
A0inv = inv(jacob(:,jdx));
A1 = jacob(:,nonzeros(DynareModel.lead_lag_incidence(1,:)));
B = jacob(:,end-number_of_shocks+1:end);

% Simulations 
for it = 2:sample_size+1
    Y(:,it) = -A0inv*(cst + A1*Y(iy1,it-1) + B*DynareOutput.exo_simul(it,:)'); 
end

DynareOutput.endo_simul = Y;