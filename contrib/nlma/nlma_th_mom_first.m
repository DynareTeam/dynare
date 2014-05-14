function moments = nlma_th_mom_first(moments,M_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nlma_th_mom_first.m
%
% This file produces the first order accurate theoretical moments, i.e.,
% the moments computed using the first order nonlinear moving average
% solution
%
% The nlma first order solution is identical to that of Dynare, so this
% file uses Dynare's solution to compute moments
%
% The mean of the first-order approximation is produced by Dynare
%
%THIS VERSION: 1.1.0 March 13, 2014
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2005-2009 Dynare Team
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

warning off MATLAB:dividebyzero

%--------------------------------------------------------------------------
% 1. Load information for moments calucation of all orders
%--------------------------------------------------------------------------
  % State variable selector
    select_state = moments.select_state;
  % Observables selector
    select_obs = moments.select_obs;
  % Number of state variables
    ns = size(select_state,2);
    
%--------------------------------------------------------------------------
% 2. Load solution from Dynare
%--------------------------------------------------------------------------
  ys  = oo_.dr.ys;        % Deterministic steady state
  ghx = oo_.dr.ghx;       % nlma's alpha
  ghu = oo_.dr.ghu;       % nlma's beta0
  Sigma_e = M_.Sigma_e;
  
%--------------------------------------------------------------------------  
% 3. Some pre-allocation
%--------------------------------------------------------------------------
  % var-cov (contemporenous covariance) matrix of the auxiliary state vector,i.e., X
    Gamma_X = cell( 1 , 1 );    
  % autocov matrix of selected variables, i.e., observables, first cell is var-cov matrix
    Gamma_obs = cell(1, options_.ar+1);
  % autocov matrix between state vector and observables, first cell is var-cov matrix
    Gamma_X_obs = cell(1, options_.ar+1);

%--------------------------------------------------------------------------
% 4. Var-cov matrices
%--------------------------------------------------------------------------
  % var-cov matrix of the state vector
    [waste, Gamma_X{1,1}] = gensylv( 1,...
                                     eye(ns),...
                                     -ghx(select_state,:),...
                                     ghx(select_state,:)',...
                                     ghu(select_state,:)*Sigma_e*ghu(select_state,:)');
  % var-cov matrix between state vector and observables
    Gamma_X_obs{1,1} = ghx(select_state,:)*Gamma_X{1,1}*ghx(select_obs,:)'...
                       + ghu(select_state,:)*Sigma_e*ghu(select_obs,:)';
  % var-cov matrix of observables
    Gamma_obs{1,1} = ghx(select_obs,:)*Gamma_X{1,1}*ghx(select_obs,:)'...
                     +  ghu(select_obs,:)*Sigma_e*ghu(select_obs,:)';    

%--------------------------------------------------------------------------
% 5. Matrix of correlation and standard deviation for observables
%--------------------------------------------------------------------------
  % Standard deviation
    moments.first_order.standard_deviation = sqrt( diag( Gamma_obs{1,1} ) );
  % Matrix of correlation
    moments.first_order.correlation = Gamma_obs{1,1}./(moments.first_order.standard_deviation*moments.first_order.standard_deviation');
   
%--------------------------------------------------------------------------                                 
% 6. Autocov matrices
%--------------------------------------------------------------------------
   if options_.ar > 0
       for j = 1: options_.ar
           % Autocov matrix of observables
             Gamma_obs{1,j+1} = ghx(select_obs,:)*Gamma_X_obs{1,j};
           % Autocov matrix between state vector and observables
             Gamma_X_obs{1,j+1} = ghx(select_state,:)*Gamma_X_obs{1,j};
           % Coefficients of autocorrelation of observables         
             moments.first_order.autocorrelation(:,j) = diag( Gamma_obs{1,j+1}./(moments.first_order.standard_deviation*moments.first_order.standard_deviation') ) ;
       end
   end
    
%--------------------------------------------------------------------------      
% 7. Mean: uses Dynare's results
%--------------------------------------------------------------------------
  moments.first_order.mean = ys(oo_.dr.order_var(select_obs),:);
  
%--------------------------------------------------------------------------
% 8. Save results for second and third order calculation
%--------------------------------------------------------------------------
  moments.first_order.Gamma_X = Gamma_X;    
  moments.first_order.Gamma_obs = Gamma_obs;
  moments.first_order.Gamma_X_obs = Gamma_X_obs;
   
    
 warning on MATLAB:dividebyzero
    
    
    
    
    
    
    
    
    
    
    
    




   




