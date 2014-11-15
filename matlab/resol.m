function [dr,info,M,options,oo] = resol(check_flag,M,options,oo)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info},@var{M},@var{options},@var{oo}] =} resol (@var{check_flag},@var{M},@var{options},@var{oo})
%! @anchor{resol}
%! @sp 1
%! Computes first and second order reduced form of the DSGE model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item check_flag
%! Integer scalar, equal to 0 if all the approximation is required, positive if only the eigenvalues are to be computed.
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item options
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @item oo
%! Matlab's structure gathering the results (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @item info==30
%! Ergodic variance can't be computed.
%! @end table
%! @sp 1
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item options
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @item oo
%! Matlab's structure gathering the results (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dynare_estimation_init}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! None.
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

if isfield(oo,'dr');
    dr = oo.dr;
end

if M.exo_nbr == 0
    oo.exo_steady_state = [] ;
end

[dr.ys,M.params,info] = evaluate_steady_state(oo.steady_state,M,options,oo,0);

if info(1)
    return
end

if options.loglinear
    % Find variables with non positive steady state.
    idx = find(dr.ys<1e-9);
    if length(idx)
        variables_with_non_positive_steady_state = M.endo_names(idx,:);
        skipline()
        fprintf('You are attempting to simulate/estimate a loglinear approximation of a model, but\n')
        fprintf('the steady state level of the following variables is not strictly positive:\n')
        for var_iter=1:length(idx)
            fprintf(' - %s (%s)\n',deblank(variables_with_non_positive_steady_state(var_iter,:)), num2str(dr.ys(idx(var_iter))))
        end
        if isestimation()
            fprintf('You should check that the priors and/or bounds over the deep parameters are such\n')
            fprintf('that the steady state levels of all the variables are strictly positive, or consider\n')
            fprintf('a linearization of the model instead of a log linearization.\n')
        else
            fprintf('You should check that the calibration of the deep parameters is such that the\n')
            fprintf('steady state levels of all the variables are strictly positive, or consider\n')
            fprintf('a linearization of the model instead of a log linearization.\n')
        end
        error('stoch_simul::resol: The loglinearization of the model cannot be performed, because the steady state is not strictly positive!')
    end
end

if options.block
    [dr,info,M,options,oo] = dr_block(dr,check_flag,M,options,oo);
else
    [dr,info] = stochastic_solvers(dr,check_flag,M,options,oo);
end
oo.dr = dr;
