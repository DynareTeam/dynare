function Ifac = mcmc_ifac(X, Nc)
% function Ifac = mcmc_ifac(X, Nc)
% Compute inefficiency factor of a MCMC sample X based on a Parzen Window
%
% INPUTS
%   X:       time series
%   Nc:      # of lags
%
% OUTPUTS
%   Ifac:       inefficiency factor of MCMC sample
%
% SPECIAL REQUIREMENTS
%   none
% ALGORITHM:
%   Inefficiency factors are computed as
%   \[
%       Ifac = 1 + 2\sum\limits_{i=1}^{Nc} {\hat \rho(i)} 
%   \]
%   where $\hat \rho(i)$ denotes the autocorrelation at lag i and the terms
%   of the sum are truncated using a Parzen window.
%   
%   For inefficiency factors, see Section 6.1 of Paolo Giordani, Michael Pitt, and Robert Kohn (2011): 
%   "Bayesian Inference for Time Series State Space Models" in : John Geweke, Gary Koop,
%   Herman van Dijk (editors): "The Oxford Handbook of Bayesian
%   Econometrics", Oxford University Press
%
%   The Parzen-Window is given by
%  \[
%   k(x) = \left\{ {\begin{array}{*{20}{c}}
%        {1 - 6{x^2} + 6|x|^3} \text{ for } 0 \leqslant |x| \leqslant \frac{1}{2}} \\
%        {2(1-|x|^3) \text{ for } \frac{1}{2} \leqslant |x| \leqslant 1} \\
%        {0 \text{ otherwise}}
%        \end{array}} \right.
%  \]
% See Donald W.K Andrews (1991): "Heteroskedasticity and autocorrelation
% consistent covariance matrix estimation", Econometrica, 59(3), p. 817-858


% Copyright (C) 2015-16 Dynare Team
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

Nc = floor(min(Nc, length(X)/2));
if mod(Nc,2),
    Nc=Nc-1;
end
AcorrXSIM = dyn_autocorr(X(:), Nc);
%
%Calculate the Parzen Weight
Parzen=zeros(Nc+1,1);
for i=1: Nc/2+1
    Parzen(i)=1 - 6*(i/Nc)^2+ 6*(i/Nc)^3;
end
for i=(Nc/2)+1: Nc+1
    Parzen(i)=2 * (1-(i/Nc))^3;
end
Parzen=Parzen';
Ifac= 1+2*sum(Parzen(:).* AcorrXSIM);
