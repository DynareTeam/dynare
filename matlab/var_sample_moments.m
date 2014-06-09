function [YtY,XtY,YtX,XtX,Y,X] = ...
    var_sample_moments(FirstObservation,LastObservation,qlag,var_trend_order,datafile,varobs,xls_sheet,xls_range)
% Computes the sample moments of a VAR model.
%
% The VAR(p) model is defined by:
%
%   y_t = \sum_{k=1}^p y_{t-k} A_k + z_t C + e_t  for t = 1,...,T  
%
% where y_t is a 1*m vector of observed endogenous variables, p is the
% number of lags, A_k is an m*m real matrix, z_t is a 1*q vector of
% exogenous (deterministic) variables, C is a q*m real matrix and
% e_t is a vector of exogenous stochastic shocks. T is the number
% of observations. The deterministic exogenous variables are assumed to 
% be a polynomial trend of order q = "var_trend_order".  
%
% We define: 
%
%  <>  Y = (y_1',y_2',...,y_T')' a T*m matrix,
%
%  <>  x_t = (y_{t-1},y_{t-2},...,y_{t-p},z_t) a 1*(mp+q) row vector, 
%
%  <>  X = (x_1',x_2',...,x_T')' a T*(mp+q) matrix, 
%
%  <>  E = (e_1',e_2',...,e_T')' a T*m matrix and
%
%  <>  A = (A_1',A_2',...,A_p',C')' an (mp+q)*m matrix of coefficients.   
%
% So that we can equivalently write the VAR(p) model using the following
% matrix representation:
%
%   Y = X * A +E
%
%
% INPUTS 
%   o FirstObservation    [integer] First observation.
%   o LastObservation     [integer] Last observation.
%   o qlag                [integer] Number of lags in the VAR model.
%   o var_trend_order     [integer] Order of the polynomial exogenous trend: 
%                                       = -1 no constant and no linear trend,
%                                       =  0 constant and no linear trend,
%                                       =  1 constant and linear trend.
%
% OUTPUTS 
%   o YtY                 [double]  Y'*Y an m*m matrix.
%   o XtY                 [double]  X'*Y an (mp+q)*m matrix. 
%   o YtX                 [double]  Y'*X an m*(mp+q) matrix.
%   o XtX                 [double]  X'*X an (mp+q)*(mp+q) matrix.
%   o Y                   [double]  Y a T*m matrix.
%   o X                   [double]  X a T*(mp+q) matrix.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2007-2009 Dynare Team
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

X = [];
Y = [];
YtY = [];
YtX = [];
XtY = [];
XtX = [];

data = read_variables(datafile,varobs,[],xls_sheet,xls_range);

if qlag >= FirstObservation
    fprintf('\nVarSampleMoments: The first used observation is %d, but there are %d lags.\n',FirstObservation,qlag)
    fprintf('VarSampleMoments: Try to increase first_obs. Given the lags, it must be at least first_obs=%d.\n',qlag+1)
    error('VarSampleMoments :: not enough data to initialize!')
    return
end

NumberOfObservations = LastObservation-FirstObservation+1;% This is T.
NumberOfVariables = size(varobs,1);% This is m.
if var_trend_order == -1% No constant no linear trend case.
    X = zeros(NumberOfObservations,NumberOfVariables*qlag);
elseif var_trend_order == 0% Constant and no linear trend case.
X = ones(NumberOfObservations,NumberOfVariables*qlag+1);
indx = NumberOfVariables*qlag+1;
elseif var_trend_order == 1;% Constant and linear trend case.
X = ones(NumberOfObservations,NumberOfVariables*qlag+2);
indx = NumberOfVariables*qlag+1:NumberOfVariables*qlag+2;
else
    error('var_sample_moments :: trend must be equal to -1,0 or 1!')
    return
end

% I build matrices Y and X  
Y = data(FirstObservation:LastObservation,:);

for t=1:NumberOfObservations
    line = t + FirstObservation-1;
    for lag = 1:qlag
        X(t,(lag-1)*NumberOfVariables+1:lag*NumberOfVariables) = data(line-lag,:);
    end
end

if (var_trend_order == 1)
    X(:,end) = transpose(1:NumberOfObservations)
end

YtY = Y'*Y;
YtX = Y'*X;
XtY = X'*Y;
XtX = X'*X;