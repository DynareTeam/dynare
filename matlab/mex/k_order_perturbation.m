% [err, g_0, g_1, g_2, g_3, derivs] = k_order_perturbation(dr,DynareModel,DynareOptions)
% computes a k_order_petrubation solution for k=1,2,3
%
% INPUTS
% dr:            struct   describing the reduced form solution of the model.
% DynareModel:   struct   jobs's parameters
% DynareOptions: struct   job's options
%
% OUTPUTS
% err:           double   err code (currently unused)
% g_0:           vector   dynare++ output. Constant effect of future volatility (in
%                         dr.order_var order) on the decision rule
%                         (=ghs2/2). Contains zero when order of
%                         approximation is 1.
% g_1:           matrix   dynare++ output. First order Taylor coefficients
%                         of decision rule. When order of approximation
%                         is 3, the. Contains both the effect of
%                         state endogenous variable and shocks. The rows
%                         are in dr.order_var order. The columns are in
%                         dr.order_var order of state endogenous
%                         variables and shocks
% g_2:           matrix   dynare++ output. Second order Taylor coefficients of decision
%                         rule. Contains both the effect of state endogenous
%                         variable and shocks. The rows are in dr.order_var
%                         order. Each row corresponds to the vectorized
%                         version of the lower triangle of the Hessian
%                         matrix. The Taylor coefficient (1/2) is
%                         included. The columns of the Hessian matrix are in
%                         dr.order_var order of state endogenous variables
%                         and shocks
% g_3:           matrix   dynare++ output. Third order Taylor coefficients of decision
%                         rule. Contains both the effect of state endogenous
%                         variable and shocks. The rows are in dr.order_var
%                         order. Each row corresponds to the vectorized
%                         version of the 3rd order derivatives tensor where each
%                         combination of variables appears only once.
%                         The Taylor coefficient (1/6) is
%                         included. Inside the tensor, the variables are in
%                         dr.order_var order of state endogenous variables
%                         and shocks 
% derivs        struct    contains the original derivatives of the
%                         decision function (ghx, ghu, ghxx, ghxu, ghuu,
%                         ghs2, ghxxx, ghxxu, ghxuu,ghuuu, ghxss, ghuss),
%                         keeping the effect of future volatility
%                         separate (in  ghs2, ghxss and ghuss). The
%                         derivatives matrices contain full versions of
%                         the Hessian matrices and 3rd order
%                         tensor. Symmetric derivatives are repeated. The
%                         Taylor coefficients (1/2 and 1/6) aren't
%                         included.
% k_order_peturbation is a compiled MEX function. It's source code is in
% dynare/mex/sources/k_order_perturbation.cc and it uses code provided by
% dynare++

% Copyright (C) 2013 Dynare Team
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
