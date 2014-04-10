function display_conditional_variance_decomposition(conditional_decomposition_array,Steps,SubsetOfVariables,M_,options_)
% This function displays the conditional variance decomposition of a given state space model
% for a subset of endogenous variables.
% 
% INPUTS 
%   conditional_decomposition_array     [matrix]   Output matrix from compute_conditional_variance_decomposition
%   Steps               [integer]     1*h vector of dates.
%   SubsetOfVariables   [integer]     1*q vector of indices.
%   M_                  [structure]   Dynare structure containing the
%                                     Model description
%   options_            [structure]   Dynare structure containing the
%                                     options
% OUTPUTS 
%   none
%
% Copyright (C) 2010-2014 Dynare Team
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

if options_.order == 2
skipline()                
disp('APPROXIMATED CONDITIONAL VARIANCE DECOMPOSITION (in percent)')
else
skipline()               
disp('CONDITIONAL VARIANCE DECOMPOSITION (in percent)')
end

vardec_i = zeros(length(SubsetOfVariables),M_.exo_nbr);

for i=1:length(Steps)
    disp(['Period ' int2str(Steps(i)) ':'])

    for j=1:M_.exo_nbr
        vardec_i(:,j) = 100*conditional_decomposition_array(:, ...
                                                          i,j);
    end
    headers = M_.exo_names;
    headers(M_.exo_names_orig_ord,:) = headers;
    headers = char(' ',headers);
    lh = size(deblank(M_.endo_names(SubsetOfVariables,:)),2)+2;
    dyntable('',headers,...
             deblank(M_.endo_names(SubsetOfVariables,:)),...
             vardec_i,lh,8,2);
end