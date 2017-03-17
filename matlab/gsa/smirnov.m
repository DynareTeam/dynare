function [H,prob,d] = smirnov(x1 , x2 , alpha, iflag )
% Smirnov test for 2 distributions
%   [H,prob,d] = smirnov(x1 , x2 , alpha, iflag )
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2012 European Commission
% Copyright (C) 2012 Dynare Team
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

if nargin<3
    alpha  =  0.05;
end
if nargin<4,
    iflag=0;
end

% empirical cdfs.
xmix= [x1;x2];
bin = [-inf ; sort(xmix) ; inf];

ncount1 = histc (x1 , bin);
ncount1 = ncount1(:);
ncount2 = histc (x2 , bin);
ncount2 = ncount2(:);

cum1  =  cumsum(ncount1)./sum(ncount1);
cum1  =  cum1(1:end-1);

cum2  =  cumsum(ncount2)./sum(ncount2);
cum2  =  cum2(1:end-1);

n1=  length(x1);
n2=  length(x2);
n =  n1*n2 /(n1+n2);

% Compute the d(n1,n2) statistics.

if iflag==0,
    d  =  max(abs(cum1 - cum2));
elseif iflag==-1
    d  =  max(cum2 - cum1);
elseif iflag==1
    d  =  max(cum1 - cum2);
end
%
% Compute P-value check H0 hypothesis
%

lam =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * d , 0);
if iflag == 0        
    j       =  [1:101]';
    prob  =  2 * sum((-1).^(j-1).*exp(-2*lam*lam*j.^2));
    
    prob=max(prob,0);
    prob=min(1,prob);
else
    prob  =  exp(-2*lam*lam);
end

H  =  (alpha >= prob);
