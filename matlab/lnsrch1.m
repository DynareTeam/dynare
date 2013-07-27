function [x,f,fvec,check]=lnsrch1(xold,fold,g,p,stpmax,func,j1,j2,varargin)
% function [x,f,fvec,check]=lnsrch1(xold,fold,g,p,stpmax,func,j1,j2,varargin)
% Computes the optimal step by minimizing the residual sum of squares
%
% INPUTS
%   xold:     actual point
%   fold:     residual sum of squares at the point xold
%   g:        gradient
%   p:        Newton direction
%   stpmax:   maximum step
%   func:     name of the function
%   j1:       equations index to be solved
%   j2:       unknowns index
%   varargin: list of arguments following j2
%
% OUTPUTS
%   x:        chosen point
%   f:        residual sum of squares value for a given x
%   fvec:     residuals vector
%   check=1:  problem of the looping which continues indefinitely
%
% 
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2010 Dynare Team
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

global options_

alf = 1e-4 ;
tolx = options_.solve_tolx;
alam = 1;

x = xold;
nn = length(j2);
summ = sqrt(sum(p.*p)) ;
if ~isfinite(summ)
    eq_number_string=[];
    for ii=1:length(j1)-1
        eq_number_string=[eq_number_string, num2str(j1(ii)), ', '];
    end
    eq_number_string=[eq_number_string, num2str(j1(end))];
    var_string=[];
    Model=evalin('base','M_');
    for ii=1:length(j2)-1
        var_string=[var_string, deblank(Model.endo_names(j2(ii),:)), ', '];
    end
    var_string=[var_string, deblank(Model.endo_names(j2(end),:))];
    fprintf('\nAn infinite element was encountered when trying to solve equation(s) %s \n',eq_number_string) 
    fprintf('with respect to the variable(s): %s.\n',var_string)
    fprintf('The values of the endogenous variables when the problem was encountered were:\n')
    for ii=1:length(xold)
        fprintf('%-s % 8.4g \n',Model.endo_names(ii,:),xold(ii));
    end
    skipline();    
    error(['Some element of Newton direction isn''t finite. Jacobian maybe' ...
           ' singular or there is a problem with initial values'])
end

if summ > stpmax
    p=p.*stpmax/summ ;
end

slope = g'*p ;

test = max(abs(p)'./max([abs(xold(j2))';ones(1,nn)])) ;
alamin = tolx/test ;

if alamin > 0.1
    alamin = 0.1;
end

while 1
    if alam < alamin
        check = 1 ;
        return
    end
    
    x(j2) = xold(j2) + (alam*p) ;
    fvec = feval(func,x,varargin{:}) ;
    fvec = fvec(j1);
    f = 0.5*fvec'*fvec ;

    if any(isnan(fvec))
        alam = alam/2 ;
        alam2 = alam ;
        f2 = f ;
        fold2 = fold ;
    else

        if f <= fold+alf*alam*slope
            check = 0;
            break ;
        else
            if alam == 1
                tmplam = -slope/(2*(f-fold-slope)) ;
            else
                rhs1 = f-fold-alam*slope ;
                rhs2 = f2-fold2-alam2*slope ;
                a = (rhs1/(alam^2)-rhs2/(alam2^2))/(alam-alam2) ;
                b = (-alam2*rhs1/(alam^2)+alam*rhs2/(alam2^2))/(alam-alam2) ;
                if a == 0
                    tmplam = -slope/(2*b) ;
                else
                    disc = (b^2)-3*a*slope ;

                    if disc < 0
                        error ('Roundoff problem in nlsearch') ;
                    else
                        tmplam = (-b+sqrt(disc))/(3*a) ;
                    end

                end

                if tmplam > 0.5*alam
                    tmplam = 0.5*alam;
                end

            end

            alam2 = alam ;
            f2 = f ;
            fold2 = fold ;
            alam = max([tmplam;(0.1*alam)]) ;
        end
    end
end

% 01/14/01 MJ lnsearch is now a separate function
% 01/12/03 MJ check for finite summ to avoid infinite loop when Jacobian
%             is singular or model is denormalized