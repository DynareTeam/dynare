function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff]  = TaRB_optimizer_wrapper(optpar,par_vector,parameterindices,TargetFun,varargin)
% function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff]  = TaRB_optimizer_wrapper(optpar,par_vector,parameterindices,TargetFun,varargin)
% Wrapper function for target function used in TaRB algorithm; reassembles
% full parameter vector before calling target function
%
% INPUTS 
%   o optpar            [double]   (p_opt*1) vector of subset of parameters to be considered 
%   o par_vector        [double]   (p*1) full vector of parameters 
%   o parameterindices  [double]   (p_opt*1) index of optpar entries in
%                                   par_vector
%   o TargetFun         [char]      string specifying the name of the objective
%                                   function (posterior kernel).
%   o varargin          [structure] other inputs of target function
% 
% OUTPUTS
%   o fval       [scalar]   value of (minus) the likelihood.
%   o info       [double]  (p*2) error code vector
%   o exit_flag  [scalar]   equal to zero if the routine return with a penalty (one otherwise).
%   o DLIK       [double]  (p*1) score vector of the likelihood.
%   o Hess       [double]  (p*p) asymptotic Hessian matrix.
%   o SteadyState [double]  Vector of doubles, steady state level for the endogenous variables.
%   o trend_coeff [double]  Matrix of doubles, coefficients of the deterministic trend in the measurement equation
% 
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

par_vector(parameterindices,:)=optpar; %reassemble parameter
[fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff] = feval(TargetFun,par_vector,varargin{:}); %call target function

