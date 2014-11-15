function [s,nu] = inverse_gamma_specification(mu,sigma,type,use_fzero_flag)
% Computes the inverse Gamma hyperparameters from the prior mean and standard deviation.

%@info:
%! @deftypefn {Function File} {[@var{s}, @var{nu} ]=} colon (@var{mu}, @var{sigma}, @var{type}, @var{use_fzero_flag})
%! @anchor{distributions/inverse_gamma_specification}
%! @sp 1
%! Computes the inverse Gamma (type 1 or 2) hyperparameters from the prior mean (@var{mu}) and standard deviation (@var{sigma}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item mu
%! Double scalar, prior mean.
%! @item sigma
%! Positive double scalar, prior standard deviation.
%! @item type
%! Integer scalar equal to one or two, type of the Inverse-Gamma distribution.
%! @item use_fzero_flag
%! Integer scalar equal to 0 (default) or 1. Use (matlab/octave's implementation of) fzero to solve for @var{nu} if equal to 1, use
%! dynare's implementation of the secant method otherwise.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item s
%! Positive double scalar (greater than two), first hypermarameter of the Inverse-Gamma prior.
%! @item nu
%! Positive double scala, second hypermarameter of the Inverse-Gamma prior.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{set_prior}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @strong{Remark:}
%! The call to the matlab's implementation of the secant method is here for testing purpose and should not be used. This routine fails
%! more often in finding an interval for nu containing a signe change because it expands the interval on both sides and eventually 
%! violates  the condition nu>2.
%! 
%! @end deftypefn
%@eod:

% Copyright (C) 2003-2012 Dynare Team
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

check_solution_flag = 1;
s = [];
nu = [];

if nargin==3
    use_fzero_flag = 0;
end

sigma2 = sigma^2;
mu2 = mu^2;

if type == 2;       % Inverse Gamma 2
    nu   = 2*(2+mu2/sigma2);
    s    = 2*mu*(1+mu2/sigma2);
elseif type == 1;   % Inverse Gamma 1
    if sigma2 < Inf
        nu = sqrt(2*(2+mu2/sigma2));
        if use_fzero_flag
            nu = fzero(@(nu)ig1fun(nu,mu2,sigma2),nu);
        else
            nu2 = 2*nu;
            nu1 = 2;
            err  = ig1fun(nu,mu2,sigma2);
            err2 = ig1fun(nu2,mu2,sigma2);
            if err2 > 0         % Too short interval.
                while nu2 < 1e12 % Shift the interval containing the root.
                    nu1  = nu2;
                    nu2  = nu2*2;
                    err2 = ig1fun(nu2,mu2,sigma2);
                    if err2<0
                        break
                    end
                end
                if err2>0
                    error('inverse_gamma_specification:: Failed in finding an interval containing a sign change! You should check that the prior variance is not too small compared to the prior mean...');
                end
            end
            % Solve for nu using the secant method.
            while abs(nu2/nu1-1) > 1e-14
                if err > 0
                    nu1 = nu;
                    if nu < nu2
                        nu = nu2;
                    else
                        nu = 2*nu;
                        nu2 = nu;
                    end
                else
                    nu2 = nu;
                end
                nu =  (nu1+nu2)/2;
                err = ig1fun(nu,mu2,sigma2);
            end
        end
        s = (sigma2+mu2)*(nu-2);
        if check_solution_flag
            if abs(log(mu)-log(sqrt(s/2))-gammaln((nu-1)/2)+gammaln(nu/2))>1e-7
                error('inverse_gamma_specification:: Failed in solving for the hyperparameters!');
            end
            if abs(sigma-sqrt(s/(nu-2)-mu^2))>1e-7
                error('inverse_gamma_specification:: Failed in solving for the hyperparameters!');
            end
        end
    else
        nu  = 2;
        s   = 2*mu2/pi;
    end
else
    error('inverse_gamma_specification: unkown type')
end

%@test:1
%$
%$ [s0,nu0] = inverse_gamma_specification(.75,.2,1,0);
%$ [s1,nu1] = inverse_gamma_specification(.75,.2,1,1);
%$ [s3,nu3] = inverse_gamma_specification(.75,.1,1,0);
%$ [s4,nu4] = inverse_gamma_specification(.75,.1,1,1);
%$ % Check the results.
%$ t(1) = dassert(s0,s1,1e-6);
%$ t(2) = dassert(nu0,nu1,1e-6);
%$ t(3) = isnan(s4);
%$ t(4) = isnan(nu4);
%$ t(5) = dassert(s3,16.240907971002265,1e-6);;
%$ t(6) = dassert(nu3,30.368398202624046,1e-6);;
%$ T = all(t);
%@eof:1