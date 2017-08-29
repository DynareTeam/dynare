function [y_out,DynareResults] =simult(y0, dr,DynareModel,DynareOptions,DynareResults)
% Simulate a DSGE model (perturbation approach).

%@info:
%! @deftypefn {Function File} {[@var{y_}, @var{DynareResults}] =} simult (@var{y0},@var{dr},@var{DynareModel},@var{DynareOptions},@var{DynareResults})
%! @anchor{simult}
%! @sp 1
%! Simulate a DSGE model (perturbation approach).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item y0
%! Vector of doubles, initial conditions.
%! @item dr
%! Matlab's structure describing decision and transition rules.
%! @item DynareModel
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_})
%! @item DynareOptions
%! Matlab's structure describing the current options (initialized by dynare, see @ref{options_}).
%! @item DynareResults
%! Matlab's structure gathering the results (see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item y_out
%! Matrix of doubles, simulated time series for all the endogenous variables (one per row).
%! @item DynareResults
%! Matlab's structure gathering the results (see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{non_linear_dsge_likelihood}, @ref{pea/pea_initialization}, @ref{stoch_simul}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{simult_}
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! If the routine is called with only one output argument, then field exo_simul (structural innovations) is not updated.
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2012 Dynare Team
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

order = DynareOptions.order;
replic = DynareOptions.simul_replic;

if replic > 1
    fname = [DynareModel.fname,'_simul'];
    fh = fopen(fname,'w+');
end

% eliminate shocks with 0 variance
i_exo_var = setdiff([1:DynareModel.exo_nbr],find(diag(DynareModel.Sigma_e) == 0));
nxs = length(i_exo_var);
DynareResults.exo_simul = zeros(DynareOptions.periods,DynareModel.exo_nbr);
chol_S = chol(DynareModel.Sigma_e(i_exo_var,i_exo_var));

for i=1:replic
    if ~isempty(DynareModel.Sigma_e)
        % we fill the shocks row wise to have the same values
        % independently of the length of the simulation
        DynareResults.exo_simul(:,i_exo_var) = randn(nxs,DynareOptions.periods)'*chol_S;
    end
    y_ = simult_(y0,dr,DynareResults.exo_simul,order);
    % elimninating initial value
    y_ = y_(:,2:end);
    if replic > 1
        fwrite(fh,y_,'float64');
    end
    if i==1
        y_out=y_;
    end
end

if replic > 1
    fclose(fh);
end