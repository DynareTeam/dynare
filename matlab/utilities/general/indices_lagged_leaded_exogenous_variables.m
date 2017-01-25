function [il,l1,ik,k1] = indices_lagged_leaded_exogenous_variables(k,M)
% returns indices of all endogenous variables split between auxiliary
% variables for lagged or leaded exogenous variables and all other ones
%
% INPUT
% k: vector of endogenous variables ID
% M: model structure
%
% OUTPUT
% il: indices of lagged or leaded variable in vector k
% l1: value of lagged or leaded variable in vector k
% ik: indices of non lagged or leaded variable in vector k
% k1: value of non lagged or leaded variable in vector k

% Copyright (C) 2011-2016 Dynare Team
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

il = [];
l1 = [];
if isempty(M.aux_vars)
    ik = 1:length(k);
    k1 = k;
else
    ik = [];
    k1 = [];
    orig_endo_nbr = M.orig_endo_nbr;
    type = [M.aux_vars.type];
    for j=1:length(k)
        if (k(j) > orig_endo_nbr)
            ty = type(k(j) - orig_endo_nbr);
            if (ty ~= 2 & ty ~= 3)
                ik = [ik; j];
                k1 = [k1; k(j)];
            else
                il = [il; j];
                l1 = [l1; k(j)];
            end
        else
            ik = [ik; j];
            k1 = [k1; k(j)];
        end
    end
end