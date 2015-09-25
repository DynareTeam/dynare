function Ifac = mcmc_ifac(X, Nc)
% function Ifac = mcmc_ifac(X, Nc)
% Compute inefficiency factor of a MCMC sample X
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

% Copyright (C) 2015 Dynare Team
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

Nc = min(Nc, length(X)/2);
AcorrXSIM = autocorr(X(:), Nc);
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
