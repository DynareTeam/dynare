function [yy, xdir, isig, lam]=log_trans_(y0,xdir0,isig,lam)

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

if nargin==4,
    % inverse transformation
    yy = (exp(y0)-lam)*isig;
    return
end

if nargin==1,
  xdir0='';
end
f=@(lam,y)gsa_skewness(log(y+lam));
isig=1;
if ~(max(y0)<0 | min(y0)>0)
  if gsa_skewness(y0)<0,
    isig=-1;
    y0=-y0;
  end
  n=hist(y0,10);
  if n(1)>20*n(end),
    try lam=fzero(f,[-min(y0)+10*eps -min(y0)+abs(median(y0))],[],y0);
    catch
      yl(1)=f(-min(y0)+10*eps,y0);
      yl(2)=f(-min(y0)+abs(median(y0)),y0);
      if abs(yl(1))<abs(yl(2))
        lam=-min(y0)+eps;
      else
        lam = -min(y0)+abs(median(y0)); %abs(100*(1+min(y0)));
      end
    end
    yy = log(y0+lam);
    xdir=[xdir0,'_logskew'];
  else
    isig=0;
    lam=0;
    yy = log(y0.^2);
    xdir=[xdir0,'_logsquared'];
  end
else
  if max(y0)<0
    isig=-1;    
    y0=-y0;
    %yy=log(-y0);
    xdir=[xdir0,'_minuslog'];
  elseif min(y0)>0
    %yy=log(y0);
    xdir=[xdir0,'_log'];
  end
  try lam=fzero(f,[-min(y0)+10*eps -min(y0)+median(y0)],[],y0);
  catch
    yl(1)=f(-min(y0)+10*eps,y0);
      yl(2)=f(-min(y0)+abs(median(y0)),y0);
    if abs(yl(1))<abs(yl(2))
      lam=-min(y0)+eps;
    else
        lam = -min(y0)+abs(median(y0)); %abs(100*(1+min(y0)));
    end
  end
  lam = max(lam,0);
  yy = log(y0+lam);
end
