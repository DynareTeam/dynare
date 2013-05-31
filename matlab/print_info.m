function print_info(info,noprint)
% Prints error messages
%
% INPUTS
%   info    [double]   vector returned by resol.m 
%   noprint [integer]  equal to 0 if the error message has to be printed. 
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2013 Dynare Team
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

if ~noprint
    switch info(1)
      case 1
        error(['The model doesn''t determine the current variables' ...
               ' uniquely'])
      case 2
        error(['The generalized Schur (QZ) decomposition failed. ' ...
               'For more information, see the documentation for Lapack function dgges: info=' ...
               int2str(info(2)) ', n=' int2str(info(3))])
      case 3
        error(['Blanchard Kahn conditions are not satisfied: no stable' ...
               ' equilibrium'])
      case 4
        error(['Blanchard Kahn conditions are not satisfied:' ...
               ' indeterminacy'])
      case 5
        error(['Blanchard Kahn conditions are not satisfied:' ...
               ' indeterminacy due to rank failure'])
      case 6
        error(['The Jacobian matrix evaluated at the steady state contains elements ' ...
               'that are not real or are infinite'])
      case 7
        error(['One of the eigenvalues is close to 0/0 (the absolute ' ...
               'value of numerator and denominator is smaller than 1e-6)'])
      case 8
        if size(info,2)>=2
          global M_;
            disp_string=deblank(M_.param_names(info(2),:));
          for ii=1:length(info)-2
            disp_string=[disp_string,', ',deblank(M_.param_names(info(2+ii),:))];
          end
          error(['The Jacobian contains NaNs because the following parameters are NaN: '...
              disp_string])
        else
          error(['The Jacobian contains NaNs'])
        end
      case 9
        error(['k_order_pert was unable to compute the solution'])           
      case 19
        error('The steadystate file did not compute the steady state')
      case 20
        error(['Impossible to find the steady state. Either the model' ...
               ' doesn''t have a steady state, there are an infinity of steady states, ' ...
               ' or the guess values are too far from the solution'])
      case 21
        error('The steady state is complex')
      case 22
        error('The steady state contains NaN or Inf')
      case 23
        error('Some updated params are complex')
      case 24
        error('Some updated params contain NaN or Inf')
      case 25
        error('The solution to the static equations is not a steady state of the dynamic model: verify that the equations tagged by [static] and [dynamic] are consistent')
      case 30
        error('Variance can''t be computed')
      case 41
        error('one (many) parameter(s) do(es) not satisfy the lower bound');
      case 42
        error('one (many) parameter(s) do(es) not satisfy the upper bound');
      case 43
        error('Covariance matrix of structural shocks is not positive definite')
      case 44 %DsgeLikelihood_hh / dsge_likelihood
        error('The covariance matrix of the measurement errors is not positive definite.');
      case 45 %DsgeLikelihood_hh / dsge_likelihood
        error('Likelihood is not a number (NaN) or a complex number');
      case 46 %DsgeLikelihood_hh / dsge_likelihood
        error('Likelihood is a complex number');
      case 47 %DsgeLikelihood_hh / dsge_likelihood
        error('Prior density is not a number (NaN)');
      case 48 %DsgeLikelihood_hh / dsge_likelihood
        error('Prior density is a complex number');
      case 49
        error('The model violates one (many) endogenous prior restriction(s)')
      case 51
        error('You are estimating a DSGE-VAR model, but the value of the dsge prior weight is too low!')
      case 52 %DsgeVarLikelihood
        error('');
      case 61 %Discretionary policy
        error(['Discretionary policy: maximum number of iterations has been reached. Procedure failed. ']);
      case 62
        error(['Discretionary policy: some eigenvalues greater than options_.qz_criterium. Model potentially unstable.']);
      case 63
        error(['Discretionary policy: NaN elements are present in the solution. Procedure failed.']);
        % Aim Code Conversions by convertAimCodeToInfo.m
      case 102
        error('Aim: roots not correctly computed by real_schur');
      case 103
        error('Aim: too many explosive roots: no stable equilibrium');
      case 135
        error('Aim: too many explosive roots, and q(:,right) is singular');
      case 104
        error('Aim: too few explosive roots: indeterminacy');
      case 145
        error('Aim: too few explosive roots, and q(:,right) is singular');
      case 105
        error('Aim: q(:,right) is singular');
      case 161
        error('Aim: too many exact shiftrights');
      case 162
        error('Aim: too many numeric shiftrights');
      case 163
        error('Aim: A is NAN or INF.')
      case 164
        error('Aim: Problem in SPEIG.')
      otherwise
        error('This case shouldn''t happen. Contact the authors of Dynare')
    end
end
