function [logged_prior_density, dlprior, d2lprior, info] = priordens(x, pshape, p6, p7, p3, p4,initialization)
% Computes a prior density for the structural parameters of DSGE models
%
% INPUTS 
%    x              [double]      vector with n elements.
%    pshape         [integer]     vector with n elements (bayestopt_.pshape).
%    p6:            [double]      vector with n elements, first  parameter of the prior distribution (bayestopt_.p6).
%    p7:            [double]      vector with n elements, second parameter of the prior distribution (bayestopt_.p7).
%    p3:            [double]      vector with n elements, lower bounds of the untruncated standard or generalized distribution
%    p4:            [double]      vector with n elements, upper bound of the untruncated standard or generalized distribution
%    initialization [integer]     if 1: initialize persistent variables
%    
% OUTPUTS 
%    logged_prior_density  [double]  scalar, log of the prior density evaluated at x.
%    info                  [double]  error code for index of Inf-prior parameter
%

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

persistent id1 id2 id3 id4 id5 id6
persistent tt1 tt2 tt3 tt4 tt5 tt6

info=0;

if nargin > 6  && initialization == 1
    % Beta indices.
    tt1 = 1;
    id1 = find(pshape==1);
    if isempty(id1)
        tt1 = 0;
    end
    % Gamma indices.
    tt2 = 1;
    id2 = find(pshape==2);
    if isempty(id2)
        tt2 = 0;
    end
    % Gaussian indices.
    tt3 = 1;
    id3 = find(pshape==3);
    if isempty(id3)
        tt3 = 0;
    end
    % Inverse-Gamma-1 indices.
    tt4 = 1;
    id4 = find(pshape==4);
    if isempty(id4)
        tt4 = 0;
    end
    % Uniform indices.
    tt5 = 1;
    id5 = find(pshape==5);
    if isempty(id5)
        tt5 = 0;
    end
    % Inverse-Gamma-2 indices.
    tt6 = 1;
    id6 = find(pshape==6);
    if isempty(id6)
        tt6 = 0;
    end
    pflag = 1;
end

logged_prior_density = 0.0;
dlprior = 0.0;
d2lprior = 0.0;

if tt1
    logged_prior_density = logged_prior_density + sum(lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1))) ;
    if isinf(logged_prior_density)
        if nargout ==4 
            info=id1(isinf(lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1))));
        end
        return
    end
    if nargout == 2,
        [tmp, dlprior(id1)]=lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1));
    elseif nargout == 3
        [tmp, dlprior(id1), d2lprior(id1)]=lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1));
    end
end

if tt2
    logged_prior_density = logged_prior_density + sum(lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2))) ;
    if isinf(logged_prior_density)
        if nargout ==4 
            info=id2(isinf(lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2))));
        end
        return
    end
    if nargout == 2,
        [tmp, dlprior(id2)]=lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2));
    elseif nargout == 3
        [tmp, dlprior(id2), d2lprior(id2)]=lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2));
    end
end

if tt3
    logged_prior_density = logged_prior_density + sum(lpdfnorm(x(id3),p6(id3),p7(id3))) ;
    if nargout == 2,
        [tmp, dlprior(id3)]=lpdfnorm(x(id3),p6(id3),p7(id3));
    elseif nargout == 3
        [tmp, dlprior(id3), d2lprior(id3)]=lpdfnorm(x(id3),p6(id3),p7(id3));
    end
end

if tt4
    logged_prior_density = logged_prior_density + sum(lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4))) ;
    if isinf(logged_prior_density)
        if nargout ==4 
            info=id4(isinf(lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4))));
        end
        return
    end
    if nargout == 2,
        [tmp, dlprior(id4)]=lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4));
    elseif nargout == 3
        [tmp, dlprior(id4), d2lprior(id4)]=lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4));
    end
end

if tt5
    if any(x(id5)-p3(id5)<0) || any(x(id5)-p4(id5)>0)
        logged_prior_density = -Inf ;
        if nargout ==4 
            info=id5((x(id5)-p3(id5)<0) || (x(id5)-p4(id5)>0));
        end
        return
    end
    logged_prior_density = logged_prior_density + sum(log(1./(p4(id5)-p3(id5)))) ;
    if nargout >1,
        dlprior(id5)=zeros(length(id5),1);
    end
    if nargout == 3
        d2lprior(id5)=zeros(length(id5),1);
    end
end

if tt6
    logged_prior_density = logged_prior_density + sum(lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6))) ;
    if isinf(logged_prior_density)
        if nargout ==4 
            info=id6(isinf(lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6))));
        end
        return
    end
    if nargout == 2,
        [tmp, dlprior(id6)]=lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6));
    elseif nargout == 3
        [tmp, dlprior(id6), d2lprior(id6)]=lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6));
    end
end

if nargout==3,
    d2lprior = diag(d2lprior);
end