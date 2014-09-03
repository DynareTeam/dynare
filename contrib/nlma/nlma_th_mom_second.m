function moments = nlma_th_mom_second(moments,M_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nlma_th_mom_second.m
%
% This file produces the second order accurate theoretical moments, i.e.,
% the moments computed using the second order nonlinear moving average
% solution
%
% The nlma second order solution is identical to that of Dynare, EXCEPT the
% constant risk correction term.
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
% 1. Compute first order accurate moments
%--------------------------------------------------------------------------
  moments = nlma_th_mom_first(moments,M_,oo_,options_);
  
%--------------------------------------------------------------------------
% 2. Load information for moments calucation of all orders
%--------------------------------------------------------------------------
  % Number of different type of variables
    nstatic = moments.nstatic;
    npred   = moments.npred;
  % State variable selector
    select_state = moments.select_state;
  % Observables selector
    select_obs = moments.select_obs;
  % Number of state variables
    ns = moments.ns;
  % Number of shocks
    ne = moments.ne;
    
%--------------------------------------------------------------------------
% 3. Load solution from Dynare and first order accurate moments
%--------------------------------------------------------------------------
  % Load solution from Dynare
    ghx  = sparse(oo_.dr.ghx);    % nlma's alpha
    ghu  = sparse(oo_.dr.ghu);    % nlma's beta0
    ghxx = sparse(oo_.dr.ghxx);   % nlma's beta22
    ghuu = sparse(oo_.dr.ghuu);   % nlma's beta00
    ghxu = sparse(oo_.dr.ghxu);   % nlma's beta20
    ghs2 = oo_.dr.ghs2;
    Sigma_e = sparse(M_.Sigma_e);
  
  % Compute nlma's y_sigma^2
    ghs2_state_nlma = (eye(npred)-ghx(nstatic+1:nstatic+npred,:))\(ghs2(nstatic+1:nstatic+npred,:));
    ghs2_nlma = [ ghx(1:nstatic,:)*ghs2_state_nlma+ghs2(1:nstatic,:)
                  ghs2_state_nlma
                  ghx(nstatic+npred+1:M_.endo_nbr,:)*ghs2_state_nlma+ghs2(nstatic+npred+1:M_.endo_nbr,:)];
    
  % Load information from first order accurate moments  
    Gamma_X_first{1,1} = sparse( moments.first_order.Gamma_X{1,1} );
    Gamma_obs_first    = moments.first_order.Gamma_obs;
    
%--------------------------------------------------------------------------  
% 4. Some presets for second order accurate moments
%--------------------------------------------------------------------------    
  % var-cov (contemporenous covariance) matrix of the auxiliary state vector,i.e., X
    Gamma_X = cell( 2 , 2 );    
  % autocov matrices of selected second order increments, i.e., observables in dy_t(2) form, first cell is var-cov matrix
    Gamma_obs = cell(1, options_.ar+1);
  % autocov matrices between state vector and selected second order increments, i.e., between X and dy_t(2), first cell is var-cov matrix
    Gamma_X_obs = cell(1, options_.ar+1);
  % autocov matrices of selected model variables, i.e.,observables in y_t(2) form, first cell is var-cov matrix
    Gamma_y_obs = cell(1, options_.ar+1);
    
%--------------------------------------------------------------------------
% 5. Build coefficient matrices
%--------------------------------------------------------------------------    
  % Theta_X
    Theta_X_11 = ghx(select_state,:);
    Theta_X_12 = 0.5*ghxx(select_state,:);
    Theta_X_22 = alt_kron(ghx(select_state,:),ghx(select_state,:));
    
  % Phi_X
    Phi_X_11 = 0.5*ghuu(select_state,:);
    Phi_X_12 = ghxu(select_state,:);
    Phi_X_21 = alt_kron(ghu(select_state,:),ghu(select_state,:));    
    temp     = alt_kron(ghx(select_state,:),ghu(select_state,:));
    Phi_X_22 = (commutation_sparse(ns,ns)+speye(ns^2))*temp;
    
  % E(Xi_t*Xi_t')
    XiXi_11  = (speye(ne^2)+commutation_sparse(ne,ne))*alt_kron(Sigma_e,Sigma_e); % Tracy and Sultan (1991), pp.344
    XiXi_22  = alt_kron(Gamma_X_first{1,1},Sigma_e);
    
  % Theta
    Theta_obs_11 = ghx(select_obs,:);
    Theta_obs_12 = 0.5*ghxx(select_obs,:);
    
  % Phi  
    Phi_obs_11 = 0.5*ghuu(select_obs,:);
    Phi_obs_12 = ghxu(select_obs,:);    
    
%--------------------------------------------------------------------------
% 2. Var-cov matrix of the auxiliary state vector, i.e., Gamma_0(X)
%--------------------------------------------------------------------------
% block 22
    temp = Phi_X_21*XiXi_11*Phi_X_21'+ Phi_X_22*XiXi_22*Phi_X_22';    
    [waste, Gamma_X{2,2}] = gensylv( 2,...
                                     eye(ns^2),...
                                     - full(Theta_X_22),...
                                     full(ghx(select_state,:))',...
                                     full(temp));
% block 21 (block 12 is equal to the transposition of block 21)
    temp = Phi_X_21*XiXi_11*Phi_X_11'+ Phi_X_22*XiXi_22*Phi_X_12';
    temp = full(Theta_X_22)*Gamma_X{2,2}*full(Theta_X_12)' + full(temp);
    [waste, Gamma_X{2,1}] = gensylv( 1,...
                                     eye(ns^2),...
                                     - full(Theta_X_22),...
                                     full(ghx(select_state,:))',...
                                     full(temp) ); 
% block 11
    temp = full(Theta_X_12)*Gamma_X{2,1}*full(Theta_X_11)';
    temp = temp + temp' + full(Theta_X_12)*Gamma_X{2,2}*full(Theta_X_12)';
    temp = Phi_X_11*XiXi_11*Phi_X_11'+ Phi_X_12*XiXi_22*Phi_X_12' + temp;
    [waste, Gamma_X{1,1}] = gensylv( 1,...
                                     eye(ns),...
                                     full(-ghx(select_state,:)),...
                                     full(ghx(select_state,:))',...
                                     full(temp));
                                                
%--------------------------------------------------------------------------
% 3. Var-cov matrix of selected model variables, i.e., Gamma_0(y)
%--------------------------------------------------------------------------
% var-cov matrix of selected second order increments, i.e., Gamma_0
   temp =  Phi_obs_11*XiXi_11*Phi_obs_11'+ Phi_obs_12*XiXi_22*Phi_obs_12';
   Gamma_obs{1,1} = ( full(Theta_obs_11)*Gamma_X{1,1}+full(Theta_obs_12)*Gamma_X{2,1} )*full(Theta_obs_11)'...
                    + ( full(Theta_obs_11)*Gamma_X{2,1}'+full(Theta_obs_12)*Gamma_X{2,2} )*full(Theta_obs_12)'...
                    + full(temp);
              
% var-cov matrix of selected model variables, i.e., Gamma_0(y)
   Gamma_y_obs{1,1} = Gamma_obs{1,1} + Gamma_obs_first{1,1};
   
%--------------------------------------------------------------------------
% 4. Var-cov matrix between state vector and selected second order 
%    increments, i.e., Gamma_0(X,dy)
%--------------------------------------------------------------------------
 temp = [ Phi_X_11*XiXi_11*Phi_obs_11'+Phi_X_12*XiXi_22*Phi_obs_12'
          Phi_X_21*XiXi_11*Phi_obs_11'+Phi_X_22*XiXi_22*Phi_obs_12'];
 Gamma_X_obs{1,1} = [ ( full(Theta_X_11)*Gamma_X{1,1}+full(Theta_X_12)*Gamma_X{2,1})*full(Theta_obs_11)'+( full(Theta_X_11)*Gamma_X{2,1}'+full(Theta_X_12)*Gamma_X{2,2})*full(Theta_obs_12)'
                       full(Theta_X_22)*Gamma_X{2,1}*full(Theta_obs_11)'+full(Theta_X_22)*Gamma_X{2,2}*full(Theta_obs_12)']...
                       + full(temp);                                     

%--------------------------------------------------------------------------
% 5. Matrix of correlation and standard deviation of selected model
%    variables
%--------------------------------------------------------------------------     
  % Standard deviation    
    moments.second_order.standard_deviation = real( sqrt( diag( Gamma_y_obs{1,1} ) ) );
  % Matrix of correlation
    moments.second_order.correlation = real( Gamma_y_obs{1,1}./(moments.second_order.standard_deviation*moments.second_order.standard_deviation') );

%--------------------------------------------------------------------------
% 6. Autocov matrices
%--------------------------------------------------------------------------
   if options_.ar > 0
        for j = 1: options_.ar
            % autocov of selected second order increments, i.e., Gamma_j
              Gamma_obs{1,j+1} = [Theta_obs_11 Theta_obs_12]*Gamma_X_obs{1,j};
            % autocov between state vector and selected second order increments, i.e., Gamma_j(X,dy)
              Gamma_X_obs{1,j+1} =  [ Theta_X_11                  Theta_X_12 
                                      sparse(ns^2,ns)             Theta_X_22 ]...
                                    *Gamma_X_obs{1,j};
            % autocov of selected model variables, i.e., Gamma_j(y)
              Gamma_y_obs{1,j+1} = Gamma_obs{1,j+1} + Gamma_obs_first{1,j+1};
            % coefficients of autocorrelation of seleced model variables         
              moments.second_order.autocorrelation(:,j) = real( diag( Gamma_y_obs{1,j+1}./(moments.second_order.standard_deviation*moments.second_order.standard_deviation') ) );
        end
   end
   
%--------------------------------------------------------------------------
% 7. Mean of selected model variables and its decomposition
%--------------------------------------------------------------------------
  % Mean
    mean_dy1state_kron2 = real( (speye(ns^2)-Theta_X_22)\(alt_kron(ghu(select_state,:),ghu(select_state,:))*vec(M_.Sigma_e)) );
    mean_dy2state       = real( (speye(ns)-ghx(select_state,:))\(0.5*ghxx(select_state,:)*mean_dy1state_kron2+0.5*ghuu(select_state,:)*vec(M_.Sigma_e)) );
    mean_dy2            = real( ghx(select_obs,:)*mean_dy2state + 0.5*ghxx(select_obs,:)*mean_dy1state_kron2 + 0.5*ghuu(select_obs,:)*vec(M_.Sigma_e) );
    moments.second_order.mean = real( full( mean_dy2 + 0.5*ghs2_nlma(select_obs,:) + moments.first_order.mean ) );
  
  % Mean decomposition: for both 2nd and 3rd accurate mean under normality,
  % Note:
  % Second order stochastic steady state = Deterministic Steady State + Risk Adjustment of Variance of Future Shock 
  %                                            Mean(2nd and 3rd Order)    Deterministic Steady State   Risk Adjustment of Variance of Future Shock      Risk Adjustment of Variance of Past Shock
    moments.second_order.mean_decomp = full( [ moments.second_order.mean  moments.first_order.mean     0.5*ghs2_nlma(select_obs,:)                      mean_dy2               ]);
    
%--------------------------------------------------------------------------
% 8. Save some terms for third order calculation
%--------------------------------------------------------------------------
  % Save nlma's y_sigma^2 
    moments.second_order.ghs2_nlma = ghs2_nlma;
  
  % Some coefficients will be recycled in third order calculation
    moments.second_order.Gamma_y_obs = Gamma_y_obs;
    
    if options_.order == 3
        moments.second_order.ghx_state_kron2     = Theta_X_22;
        moments.second_order.ghx_state_kron3     = alt_kron( Theta_X_22 , ghx(select_state,:) );
        moments.second_order.ghu_state_kron2     = Phi_X_21;
        moments.second_order.ghu_state_kron3     = alt_kron(ghu(select_state,:),Phi_X_21);
        moments.second_order.ee_kron3            = alt_kron(alt_kron(Sigma_e,Sigma_e),Sigma_e);
        moments.second_order.mean_dy1state_kron2 = mean_dy1state_kron2;
        moments.second_order.mean_dy2state       = mean_dy2state;
        moments.second_order.mean_dy2            = mean_dy2;
        moments.second_order.Gamma_X             = Gamma_X;
        moments.second_order.Gamma_obs{1,1}      = Gamma_obs{1,1};        
    end


 warning on MATLAB:dividebyzero




   




