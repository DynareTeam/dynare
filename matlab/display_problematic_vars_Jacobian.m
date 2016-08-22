function []=display_problematic_vars_Jacobian(problemrow,problemcol,M_,x,type,caller_string)
% []=display_problematic_vars_Jacobian(problemrow,problemcol,M_,ys,caller_string)
% print the equation numbers and variables associated with problematic entries 
% of the Jacobian 
%
% INPUTS
%   problemrow      [vector] rows associated with problematic entries
%   problemcol      [vector] columns associated with problematic entries
%   M_              [matlab structure] Definition of the model.           
%   x               [vector] point at which the Jacobian was evaluated
%   type            [string] 'static' or 'dynamic' depending on the type of
%                               Jacobian
%   caller_string   [string] contains name of calling function for printing 
%    
% OUTPUTS
%   none.
%  

% Copyright (C) 2014-16 Dynare Team
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

skipline();
if nargin<6
    caller_string='';
end
aux_eq_nbr=M_.eq_nbr-M_.orig_eq_nbr;
if strcmp(type,'dynamic')
    for ii=1:length(problemrow)
        if problemcol(ii)>max(M_.lead_lag_incidence)
            var_row=2;
            var_index=problemcol(ii)-max(max(M_.lead_lag_incidence));
        else
            [var_row,var_index]=find(M_.lead_lag_incidence==problemcol(ii));
        end
        if var_row==2
            type_string='';
        elseif var_row==1
            type_string='lag of';
        elseif var_row==3;
            type_string='lead of';
        end
        if problemcol(ii)<=max(max(M_.lead_lag_incidence)) && var_index<=M_.orig_endo_nbr
            if problemrow(ii)<=aux_eq_nbr
                eq_nbr=problemrow(ii);
                fprintf('Derivative of Auxiliary Equation %d with respect to %s Variable %s  (initial value of %s: %g) \n',eq_nbr,type_string,deblank(M_.endo_names(var_index,:)),deblank(M_.endo_names(var_index,:)),x(var_index))
            else
                eq_nbr=problemrow(ii)-aux_eq_nbr;
                fprintf('Derivative of Equation %d with respect to %s Variable %s  (initial value of %s: %g) \n',eq_nbr,type_string,deblank(M_.endo_names(var_index,:)),deblank(M_.endo_names(var_index,:)),x(var_index))
            end
        elseif problemcol(ii)<=max(max(M_.lead_lag_incidence)) && var_index>M_.orig_endo_nbr %auxiliary vars
            if M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).type ==6 %Ramsey Lagrange Multiplier 
                if problemrow(ii)<=aux_eq_nbr
                    eq_nbr=problemrow(ii);
                    fprintf('Derivative of Auxiliary Equation %d with respect to %s of Langrange multiplier of equation %s (initial value: %g) \n',eq_nbr,type_string,M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).eq_nbr,x(problemcol(ii)))            
                else
                    eq_nbr=problemrow(ii)-aux_eq_nbr;
                    fprintf('Derivative of Equation %d with respect to %s of Langrange multiplier of equation %s (initial value: %g) \n',eq_nbr,type_string,M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).eq_nbr,x(problemcol(ii)))            
                end
            else
                if problemrow(ii)<=aux_eq_nbr
                    eq_nbr=problemrow(ii);
                    orig_var_index=M_.aux_vars(1,var_index-M_.orig_endo_nbr).orig_index;
                    fprintf('Derivative of Auxiliary Equation %d with respect to %s Variable %s  (initial value of %s: %g) \n',eq_nbr,type_string,deblank(M_.endo_names(orig_var_index,:)),deblank(M_.endo_names(orig_var_index,:)),x(orig_var_index))            
                else
                    eq_nbr=problemrow(ii)-aux_eq_nbr;
                    orig_var_index=M_.aux_vars(1,var_index-M_.orig_endo_nbr).orig_index;
                    fprintf('Derivative of Equation %d with respect to %s Variable %s  (initial value of %s: %g) \n',eq_nbr,type_string,deblank(M_.endo_names(orig_var_index,:)),deblank(M_.endo_names(orig_var_index,:)),x(orig_var_index))            
                end
            end
        elseif problemcol(ii)>max(max(M_.lead_lag_incidence)) && var_index<=M_.exo_nbr
            if problemrow(ii)<=aux_eq_nbr
                eq_nbr=problemrow(ii);
                fprintf('Derivative of Auxiliary Equation %d with respect to %s shock %s \n',eq_nbr,type_string,deblank(M_.exo_names(var_index,:)));
            else
                eq_nbr=problemrow(ii)-aux_eq_nbr;
                fprintf('Derivative of Equation %d with respect to %s shock %s \n',eq_nbr,type_string,deblank(M_.exo_names(var_index,:)));
            end            
        else
            error('display_problematic_vars_Jacobian:: The error should not happen. Please contact the developers')
        end
    end
    fprintf('\n%s  The problem most often occurs, because a variable with\n',caller_string)
    fprintf('%s  exponent smaller than 1 has been initialized to 0. Taking the derivative\n',caller_string)
    fprintf('%s  and evaluating it at the steady state then results in a division by 0.\n',caller_string)
    fprintf('%s  If you are using model-local variables (# operator), check their values as well.\n',caller_string)
elseif strcmp(type,'static')
    for ii=1:length(problemrow)
        if problemcol(ii)<=M_.orig_endo_nbr
            if problemrow(ii)<=aux_eq_nbr
                eq_nbr=problemrow(ii);
                fprintf('Derivative of Auxiliary Equation %d with respect to Variable %s  (initial value of %s: %g) \n',eq_nbr,deblank(M_.endo_names(problemcol(ii),:)),deblank(M_.endo_names(problemcol(ii),:)),x(problemcol(ii)))
            else
                eq_nbr=problemrow(ii)-aux_eq_nbr;
                fprintf('Derivative of Equation %d with respect to Variable %s  (initial value of %s: %g) \n',eq_nbr,deblank(M_.endo_names(problemcol(ii),:)),deblank(M_.endo_names(problemcol(ii),:)),x(problemcol(ii)))
            end
        else %auxiliary vars
            if M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).type ==6 %Ramsey Lagrange Multiplier 
                if problemrow(ii)<=aux_eq_nbr
                    eq_nbr=problemrow(ii);
                    fprintf('Derivative of Auxiliary Equation %d with respect to Lagrange multiplier of equation %d (initial value: %g) \n',eq_nbr,M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).eq_nbr,x(problemcol(ii)))            
                else
                    eq_nbr=problemrow(ii)-aux_eq_nbr;
                    fprintf('Derivative of Equation %d with respect to Lagrange multiplier of equation %d (initial value: %g) \n',eq_nbr,M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).eq_nbr,x(problemcol(ii)))            
                end
            else
            if problemrow(ii)<=aux_eq_nbr
                eq_nbr=problemrow(ii);
                orig_var_index=M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).orig_index;
                fprintf('Derivative of Auxiliary Equation %d with respect to Variable %s  (initial value of %s: %g) \n',eq_nbr,deblank(M_.endo_names(orig_var_index,:)),deblank(M_.endo_names(orig_var_index,:)),x(problemcol(ii)))            
            else
                eq_nbr=problemrow(ii)-aux_eq_nbr;
                orig_var_index=M_.aux_vars(1,problemcol(ii)-M_.orig_endo_nbr).orig_index;
                fprintf('Derivative of Equation %d with respect to Variable %s  (initial value of %s: %g) \n',eq_nbr,deblank(M_.endo_names(orig_var_index,:)),deblank(M_.endo_names(orig_var_index,:)),x(problemcol(ii)))            
            end
            end
        end
    end
    fprintf('\n%s  The problem most often occurs, because a variable with\n',caller_string)
    fprintf('%s  exponent smaller than 1 has been initialized to 0. Taking the derivative\n',caller_string)
    fprintf('%s  and evaluating it at the steady state then results in a division by 0.\n',caller_string)
    fprintf('%s  If you are using model-local variables (# operator), check their values as well.\n',caller_string)
else
    error('Unknown Type')    
end