function moments = nlma_th_mom_third(moments,M_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nlma_th_mom_third.m
%
% This file produces the third order accurate theoretical moments, i.e.,
% the moments computed using the third order nonlinear moving average
% solution, and a variance decomposition
%
% The nlma third order solution is identical to that of Dynare, EXCEPT the
% constant and time varying risk correction terms.
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
% 1. Compute first and second order accurate moments
%--------------------------------------------------------------------------
  moments = nlma_th_mom_second(moments,M_,oo_,options_);
  
%--------------------------------------------------------------------------
% 2. Load information for all three orders
%--------------------------------------------------------------------------
  % Number of different type of variables
    npred = moments.npred;
  % State variable selector
    select_state = moments.select_state;
  % Observables selector
    select_obs = moments.select_obs;
  % Number of state variables
    ns = moments.ns;
  % Number of observables
    nobs = moments.nobs;
  % Number of shocks
    ne = moments.ne;
    
%--------------------------------------------------------------------------
% 3. Load solution from Dynare and first and second order accurate moments
%--------------------------------------------------------------------------
  % Load solution from Dynare
    ghx   = sparse(oo_.dr.ghx);    % nlma's alpha
    ghu   = sparse(oo_.dr.ghu);    % nlma's beta0
    ghxx  = sparse(oo_.dr.ghxx);   % nlma's beta22
    ghuu  = sparse(oo_.dr.ghuu);   % nlma's beta00
    ghxu  = sparse(oo_.dr.ghxu);   % nlma's beta20    
    ghxxx = sparse(oo_.dr.ghxxx);  % nlma's beta333_1
    ghxxu = sparse(oo_.dr.ghxxu);  % nlma's beta330_1
    ghxuu = sparse(oo_.dr.ghxuu);  % nlma's beta300
    ghuuu = sparse(oo_.dr.ghuuu);  % nlma's beta000
    ghxss = sparse(oo_.dr.ghxss);
    ghuss = sparse(oo_.dr.ghuss);
    Sigma_e = sparse(M_.Sigma_e);
    
  % Compute some nlma coefficients
    % nlma's y_sigma^2, already known from second order accurate moments
    % caculation
      ghs2_nlma = sparse( moments.second_order.ghs2_nlma );       
    % nlma's y_sigma^2e, or beta_sigma2_0          
      ghuss_nlma = ghuss + ghxu*alt_kron(ghs2_nlma(select_state,:),eye(ne));
    % nlma's y_sigma^2y^state, or beta_sigma2_1
      ghxss_nlma = ghxss + ghxx*alt_kron(ghs2_nlma(select_state,:),eye(npred));
  
  % Load information from first and second order accurate moments  
    Gamma_X_first{1,1}  = sparse( moments.first_order.Gamma_X{1,1} );
    Gamma_X_second      = moments.second_order.Gamma_X;
    Gamma_obs_first{1,1}  = moments.first_order.Gamma_obs{1,1};
    Gamma_obs_second{1,1} = moments.second_order.Gamma_obs{1,1};
    Gamma_y_obs_second  = moments.second_order.Gamma_y_obs;    
    mean_dy1state_kron2 = moments.second_order.mean_dy1state_kron2;
    mean_dy2state       = moments.second_order.mean_dy2state;
    ghx_state_kron2     = sparse( moments.second_order.ghx_state_kron2 );
    ghx_state_kron3     = sparse( moments.second_order.ghx_state_kron3 );
    ghu_state_kron2     = sparse( moments.second_order.ghu_state_kron2 );
    ghu_state_kron3     = sparse( moments.second_order.ghu_state_kron3 );
    ee_kron3            = sparse( moments.second_order.ee_kron3 );

%--------------------------------------------------------------------------  
% 4. Some presets for third order accurate moments
%--------------------------------------------------------------------------
  % var-cov (contemporaneous covariance) matrix of the auxiliary state vector, i.e., X
    Gamma_X = cell( 4 , 4  );
  % autocov matrices of selected third order increments, i.e., observables in dy_t(3) form, first cel is var-cov matrix
    Gamma_obs = cell(1, options_.ar+1);
  % autocovariance matrices between state vector and selected third order increments, i.e., between X and dy_t(3), first cell is var-cov matrix
    Gamma_X_obs = cell(4, options_.ar+1);
  % autocovariance matrices of selected model variables, i.e., observables in y_t(3) form, first cell is var-cov matrix
    Gamma_y_obs = cell(1, options_.ar+1); 
  % autocovariance matrices between selected third and first order increments, i.e., between dy_t(1) and dy_t(3), first cell is var-cov matrix
    Gamma_13_obs = cell(1, options_.ar+1); 
    
%--------------------------------------------------------------------------
% 1. Build coefficient matrices
%--------------------------------------------------------------------------   
  % Phi_X      
    Phi_X_11 = (1/6)*ghuuu(select_state,:); 
    Phi_X_12 = 0.5*ghxxu(select_state,:);
    Phi_X_13 = 0.5*ghxuu(select_state,:);
    Phi_X_14 = ghxu(select_state,:);
    Phi_X_15 = 0.5*ghuss_nlma(select_state,:);
    Phi_X_21 = alt_kron( 0.5*ghuu(select_state,:) , ghu(select_state,:)  );    
    Phi_X_22 = alt_kron( ghxu(select_state,:) , ghx(select_state,:) )*commutation_sparse(ns*ne,ns) + alt_kron(0.5*ghxx(select_state,:),ghu(select_state,:));
    Phi_X_23 = alt_kron( 0.5*ghuu(select_state,:) , ghx(select_state,:) )*commutation_sparse(ne^2,ns) + alt_kron(ghxu(select_state,:),ghu(select_state,:) );
    Phi_X_24 = alt_kron( ghx(select_state,:) , ghu(select_state,:) ); 
    Phi_X_31 = ghu_state_kron3;   
    temp = (alt_kron(commutation_sparse(ns,ns),speye(ns))+speye(ns^3))*commutation_sparse(ns^2,ns)+speye(ns^3);
    Phi_X_32 = temp*alt_kron(ghx_state_kron2 , ghu(select_state,:));
    temp = commutation_sparse(ns^2,ns)+alt_kron(commutation_sparse(ns,ns),speye(ns))+speye(ns^3);
    Phi_X_33 = temp*alt_kron(ghx(select_state,:) , ghu_state_kron2);
    Phi_X_45 = ghu(select_state,:);
    
  % Theta_X
    Theta_X_11 = ghx(select_state,:);
    Theta_X_12 = ghxx(select_state,:);
    Theta_X_13 = (1/6)*ghxxx(select_state,:);
    Theta_X_14 = 0.5*ghxuu(select_state,:)*alt_kron( speye(ns),vec(Sigma_e) )+ 0.5*ghxss_nlma(select_state,:);
    Theta_X_22 = ghx_state_kron2;
    Theta_X_23 = alt_kron( 0.5*ghxx(select_state,:) , ghx(select_state,:) );
    Theta_X_24 = Phi_X_23*alt_kron( speye(ns),vec(Sigma_e) ); 
    Theta_X_33 = ghx_state_kron3;    
    Theta_X_34 = Phi_X_33*alt_kron( speye(ns),vec(Sigma_e) );
    Theta_X_44 = ghx(select_state,:); 
    
  % Theta_obs
    Theta_obs_11 = ghx(select_obs,:);
    Theta_obs_12 = ghxx(select_obs,:);
    Theta_obs_13 = (1/6)*ghxxx(select_obs,:);
    Theta_obs_14 = 0.5*ghxuu(select_obs,:)*alt_kron(speye(ns),vec(Sigma_e)) + 0.5*ghxss_nlma(select_obs,:);    
    
  % Phi_obs
    Phi_obs_11 = (1/6)*ghuuu(select_obs,:); 
    Phi_obs_12 = 0.5*ghxxu(select_obs,:);
    Phi_obs_13 = 0.5*ghxuu(select_obs,:);
    Phi_obs_14 = ghxu(select_obs,:);        
    Phi_obs_15 = 0.5*ghuss_nlma(select_obs,:);
    
  % E(Xi*Xi'), it is symmetric     
    XiXi_55 = Sigma_e;
    XiXi_54 = alt_kron(mean_dy2state',Sigma_e);
    XiXi_52 = alt_kron(mean_dy1state_kron2',Sigma_e);    
    temp = alt_kron(vec(Sigma_e),Sigma_e);
    XiXi_51 = ( alt_kron(Sigma_e,vec(Sigma_e))+temp+alt_kron(speye(ne),commutation_sparse(ne,ne))*temp )'; % Jinadasa and Tracy (1986, p.404)
    temp = Gamma_X_second{1,1} + mean_dy2state*mean_dy2state';
    XiXi_44 = alt_kron( temp ,Sigma_e);
    temp = Gamma_X_second{2,1}' + mean_dy2state*mean_dy1state_kron2';
    XiXi_42 = alt_kron( temp , Sigma_e );    
    XiXi_41 = alt_kron(mean_dy2state,XiXi_51);    
    temp = (speye(ne^2)+commutation_sparse(ne,ne))*alt_kron(Sigma_e,Sigma_e); % Tracy and Sultan (1993, p.344)
    XiXi_33 = alt_kron(  Gamma_X_first{1,1}, temp );   
    temp = Gamma_X_second{2,2}+mean_dy1state_kron2*mean_dy1state_kron2';
    XiXi_22 = alt_kron( temp ,Sigma_e  );    
    XiXi_21 = alt_kron(  mean_dy1state_kron2,XiXi_51  );    
    K_for_six_mom = speye(ne^3) + commutation_sparse(ne^2,ne) + commutation_sparse(ne,ne^2); % Tracy and Sultan (1993, p.344-345)
    temp = K_for_six_mom + alt_kron(speye(ne),commutation_sparse(ne,ne)) + alt_kron(commutation_sparse(ne,ne),speye(ne)) + commutation_sparse(ne,ne^2)*alt_kron(commutation_sparse(ne,ne),speye(ne)); 
    XiXi_11 = ee_kron3*temp+K_for_six_mom*alt_kron( vec(Sigma_e)*vec(Sigma_e)',Sigma_e )*K_for_six_mom;
     
%--------------------------------------------------------------------------
% 2. Var-cov matrix of the auxiliary state vector, i.e., Gamma_0(X)
%--------------------------------------------------------------------------    
  % block 44
      Gamma_X{4,4} = full(Gamma_X_first{1,1});
  % block 43     
      temp = Phi_X_45*( XiXi_51*Phi_X_31'+XiXi_52*Phi_X_32' ); 
      temp = Theta_X_44*Gamma_X{4,4}*Theta_X_34' + temp; 
      [waste, Gamma_X{4,3}] = gensylv( 3,...
                                       eye(ns),...
                                       full( -ghx(select_state,:) ),...
                                       full( ghx(select_state,:)' ),...
                                       full( temp )  );
  % block 42           
      temp = Phi_X_45*( XiXi_51*Phi_X_21'+XiXi_52*Phi_X_22'+XiXi_54*Phi_X_24' ); 
      temp = Theta_X_44*( Gamma_X{4,3}*Theta_X_23'+Gamma_X{4,4}*Theta_X_24' )+ temp;   
      [waste, Gamma_X{4,2}] = gensylv( 2,...
                                       eye(ns),...
                                       full( -ghx(select_state,:) ),...
                                       full( ghx(select_state,:)' ),...
                                       full(temp) );      
  % block 41
      temp = Phi_X_45*( XiXi_51*Phi_X_11'+XiXi_52*Phi_X_12'+XiXi_54*Phi_X_14'+XiXi_55*Phi_X_15' );
      temp = Theta_X_44*(Gamma_X{4,2}*Theta_X_12'+Gamma_X{4,3}*Theta_X_13'+Gamma_X{4,4}*Theta_X_14')+ temp;    
      [waste, Gamma_X{4,1}] = gensylv( 1,...
                                       eye(ns),...
                                       full( -ghx(select_state,:) ),...
                                       full( ghx(select_state,:)' ),...
                                       full(temp) );      
 % block 33
      temp =  ( full(Phi_X_31)*full(XiXi_11)+full(Phi_X_32*XiXi_21) )*full(Phi_X_31')...
             +( full(Phi_X_31)*full(XiXi_21')+full(Phi_X_32)*full(XiXi_22) )*full(Phi_X_32') + full(Phi_X_33)*full(XiXi_33)*full(Phi_X_33');
      temp =   full(Theta_X_34)*Gamma_X{4,3}*full(Theta_X_33')...
             +( full(Theta_X_33)*Gamma_X{4,3}'+full(Theta_X_34)*Gamma_X{4,4} )*full(Theta_X_34') + temp;
      if ns >= 8
          % In the case the number of state variables is equal to or larger
          % than 8, disclyap_kron_3.m is faster than Dynare's gensylv.m
          Gamma_X{3,3} = disclyap_kron_3( ghx(select_state,:),temp,1 );
      else
          [waste, Gamma_X{3,3}] = gensylv( 3,...
                                           eye(ns^3),...
                                           full( -Theta_X_33 ),...
                                           full( ghx(select_state,:)' ),...
                                           temp );
      end                                                        
 % block 32
      temp =  ( full(Phi_X_31)*full(XiXi_11)+full(Phi_X_32)*full(XiXi_21) )*full(Phi_X_21') + ( full(Phi_X_31)*full(XiXi_21')+full(Phi_X_32)*full(XiXi_22) )*full(Phi_X_22')...
             + full(Phi_X_33)*full(XiXi_33)*full(Phi_X_23') + ( full(Phi_X_31)*full(XiXi_41')+full(Phi_X_32)*full(XiXi_42') )*full(Phi_X_24');
      temp =  full(Theta_X_34)*Gamma_X{4,2}*full(Theta_X_22') + ( full(Theta_X_33)*Gamma_X{3,3}+full(Theta_X_34)*Gamma_X{4,3} )*full(Theta_X_23')...
             + ( full(Theta_X_33)*Gamma_X{4,3}'+full(Theta_X_34)*Gamma_X{4,4} )*full(Theta_X_24') + temp;      
      [waste, Gamma_X{3,2}] = gensylv( 3,...
                                       eye(ns^2),...
                                       full( -Theta_X_22 ),...
                                       full( ghx(select_state,:)' ),...
                                       temp' );
      Gamma_X{3,2} = Gamma_X{3,2}';
 % block 31
      temp =  ( full(Phi_X_31)*full(XiXi_11)+full(Phi_X_32)*full(XiXi_21) )*full(Phi_X_11') + ( full(Phi_X_31)*full(XiXi_21')+full(Phi_X_32)*full(XiXi_22) )*full(Phi_X_12')...
             + full(Phi_X_33)*full(XiXi_33)*full(Phi_X_13') + ( full(Phi_X_31)*full(XiXi_41')+full(Phi_X_32)*full(XiXi_42') )*full(Phi_X_14') + ( full(Phi_X_31)*full(XiXi_51')+full(Phi_X_32)*full(XiXi_52') )*full(Phi_X_15');
      temp =  full(Theta_X_34)*Gamma_X{4,1}*full(Theta_X_11') + ( full(Theta_X_33)*Gamma_X{3,2}+full(Theta_X_34)*Gamma_X{4,2} )*full(Theta_X_12')...
             +( full(Theta_X_33)*Gamma_X{3,3}+full(Theta_X_34)*Gamma_X{4,3})*full(Theta_X_13') + ( full(Theta_X_33)*Gamma_X{4,3}'+full(Theta_X_34)*Gamma_X{4,4})*full(Theta_X_14') + temp;
      [waste, Gamma_X{3,1}] = gensylv( 3,...
                                       eye(ns),...
                                       full( -ghx(select_state,:) ),...
                                       full( ghx(select_state,:)' ),...
                                       temp' );                                                      
       Gamma_X{3,1} = Gamma_X{3,1}';
 % block 22
       temp = ( full(Phi_X_21)*full(XiXi_11)+full(Phi_X_22)*full(XiXi_21)+full(Phi_X_24)*full(XiXi_41) )*full(Phi_X_21') + ( full(Phi_X_21)*full(XiXi_21')+full(Phi_X_22)*full(XiXi_22)+full(Phi_X_24)*full(XiXi_42) )*full(Phi_X_22') + full(Phi_X_23)*full(XiXi_33)*full(Phi_X_23') + ( full(Phi_X_21)*full(XiXi_41')+full(Phi_X_22)*full(XiXi_42')+full(Phi_X_24)*full(XiXi_44) )*full(Phi_X_24'); 
       temp = ( full(Theta_X_23)*Gamma_X{3,2}+full(Theta_X_24)*Gamma_X{4,2})*full(Theta_X_22')+( full(Theta_X_22)*Gamma_X{3,2}'+full(Theta_X_23)*Gamma_X{3,3}+full(Theta_X_24)*Gamma_X{4,3})*full(Theta_X_23')...
             +( full(Theta_X_22)*Gamma_X{4,2}'+full(Theta_X_23)*Gamma_X{4,3}'+full(Theta_X_24)*Gamma_X{4,4})*full(Theta_X_24') + temp;
       [waste, Gamma_X{2,2}] = gensylv( 2,...
                                        eye(ns^2),...
                                        full( -Theta_X_22 ),...
                                        full( ghx(select_state,:)' ),...
                                        temp );  
 % block 21
       temp = (full(Phi_X_21)*full(XiXi_11)+full(Phi_X_22)*full(XiXi_21)+full(Phi_X_24)*full(XiXi_41))*full(Phi_X_11') + (full(Phi_X_21)*full(XiXi_21')+full(Phi_X_22)*full(XiXi_22)+full(Phi_X_24)*full(XiXi_42))*full(Phi_X_12') + full(Phi_X_23)*full(XiXi_33)*full(Phi_X_13') + (full(Phi_X_21)*full(XiXi_41')+full(Phi_X_22)*full(XiXi_42')+full(Phi_X_24)*full(XiXi_44))*full(Phi_X_14') + (full(Phi_X_21)*full(XiXi_51')+full(Phi_X_22)*full(XiXi_52')+full(Phi_X_24)*full(XiXi_54'))*full(Phi_X_15');
       temp =  ( full(Theta_X_23)*Gamma_X{3,1}+full(Theta_X_24)*Gamma_X{4,1} )*full(Theta_X_11')+( full(Theta_X_22)*Gamma_X{2,2} + full(Theta_X_23)*Gamma_X{3,2} + full(Theta_X_24)*Gamma_X{4,2})*full(Theta_X_12')...
              +( full(Theta_X_22)*Gamma_X{3,2}'+ full(Theta_X_23)*Gamma_X{3,3} + full(Theta_X_24)*Gamma_X{4,3})*full(Theta_X_13')+( full(Theta_X_22)*Gamma_X{4,2}'+ full(Theta_X_23)*Gamma_X{4,3}'+ full(Theta_X_24)*Gamma_X{4,4})*full(Theta_X_14') + temp; 
       [waste, Gamma_X{2,1}] = gensylv( 1,...
                                        eye(ns^2),...
                                        full( -Theta_X_22 ),...
                                        full( ghx(select_state,:)' ),...
                                        temp );
 % block 11
       temp = (full(Phi_X_11)*full(XiXi_11)+full(Phi_X_12)*full(XiXi_21)+full(Phi_X_14)*full(XiXi_41)+full(Phi_X_15)*full(XiXi_51))*full(Phi_X_11') + (full(Phi_X_11)*full(XiXi_21')+full(Phi_X_12)*full(XiXi_22)+full(Phi_X_14)*full(XiXi_42)+full(Phi_X_15)*full(XiXi_52))*full(Phi_X_12') + full(Phi_X_13)*full(XiXi_33)*full(Phi_X_13') + (full(Phi_X_11)*full(XiXi_41')+full(Phi_X_12)*full(XiXi_42')+full(Phi_X_14)*full(XiXi_44)+full(Phi_X_15)*full(XiXi_54))*full(Phi_X_14') + (full(Phi_X_11)*full(XiXi_51')+full(Phi_X_12)*full(XiXi_52')+full(Phi_X_14)*full(XiXi_54')+full(Phi_X_15)*full(XiXi_55))*full(Phi_X_15'); 
       temp = (full(Theta_X_12)*Gamma_X{2,1}+full(Theta_X_13)*Gamma_X{3,1}+full(Theta_X_14)*Gamma_X{4,1})*full(Theta_X_11') + (full(Theta_X_11)*Gamma_X{2,1}'+full(Theta_X_12)*Gamma_X{2,2}+full(Theta_X_13)*Gamma_X{3,2}+full(Theta_X_14)*Gamma_X{4,2})*full(Theta_X_12')...
             +(full(Theta_X_11)*Gamma_X{3,1}'+full(Theta_X_12)*Gamma_X{3,2}'+full(Theta_X_13)*Gamma_X{3,3}+full(Theta_X_14)*Gamma_X{4,3})*full(Theta_X_13') +(full(Theta_X_11)*Gamma_X{4,1}'+full(Theta_X_12)*Gamma_X{4,2}'+full(Theta_X_13)*Gamma_X{4,3}'+full(Theta_X_14)*Gamma_X{4,4})*full(Theta_X_14') + temp; 
       [waste, Gamma_X{1,1}] = gensylv( 1,...
                                        eye(ns),...
                                        full( -ghx(select_state,:) ),...
                                        full( ghx(select_state,:)' ),...
                                        temp );
  
%--------------------------------------------------------------------------
% 3. Var-cov matrix of selected model variables, i.e., Gamma_0(y)
%--------------------------------------------------------------------------   
% var-cov matrix of selected third order increments, i.e., Gamma_0
  temp = (full(Phi_obs_11)*full(XiXi_11) + full(Phi_obs_12)*full(XiXi_21)+full(Phi_obs_14)*full(XiXi_41)+full(Phi_obs_15)*full(XiXi_51))*full(Phi_obs_11')...
         +(full(Phi_obs_12)*full(XiXi_22)+ full(Phi_obs_14)*full(XiXi_42) + full(Phi_obs_15)*full(XiXi_52) + full(Phi_obs_11)*full(XiXi_21'))*full(Phi_obs_12')...  
         +(full(Phi_obs_14)*full(XiXi_44) + full(Phi_obs_15)*full(XiXi_54) + full(Phi_obs_11)*full(XiXi_41') + full(Phi_obs_12)*full(XiXi_42'))*full(Phi_obs_14')...
         +(full(Phi_obs_15)*full(XiXi_55) + full(Phi_obs_11)*full(XiXi_51') + full(Phi_obs_12)*full(XiXi_52') + full(Phi_obs_14)*full(XiXi_54'))*full(Phi_obs_15')+ full(Phi_obs_13)*full(XiXi_33)*full(Phi_obs_13');
  Gamma_obs{1,1} =   ( full(Theta_obs_11)*Gamma_X{1,1}+full(Theta_obs_12)*Gamma_X{2,1}+full(Theta_obs_13)*Gamma_X{3,1}+full(Theta_obs_14)*Gamma_X{4,1} )*full(Theta_obs_11') + ( full(Theta_obs_12)*Gamma_X{2,2}+full(Theta_obs_11)*Gamma_X{2,1}'+full(Theta_obs_13)*Gamma_X{3,2}+full(Theta_obs_14)*Gamma_X{4,2} )*full(Theta_obs_12')... 
                     + ( full(Theta_obs_13)*Gamma_X{3,3}+full(Theta_obs_11)*Gamma_X{3,1}'+full(Theta_obs_12)*Gamma_X{3,2}'+full(Theta_obs_14)*Gamma_X{4,3} )*full(Theta_obs_13') + ( full(Theta_obs_11)*Gamma_X{4,1}'+full(Theta_obs_14)*Gamma_X{4,4}+full(Theta_obs_12)*Gamma_X{4,2}'+full(Theta_obs_13)*Gamma_X{4,3}' )*full(Theta_obs_14') + temp;

% var-cov matrix of selected model variables, i.e., Gamma_0(y)
   % var-cov matrix between selected first and third order increments, i.e., between dy_t(1) and dy_t(3)   
      Gamma_13_obs{1,1} = full(ghx(select_obs,:))*( Gamma_X{4,1}*full(Theta_obs_11')+Gamma_X{4,2}*full(Theta_obs_12')+ Gamma_X{4,3}*full(Theta_obs_13')+Gamma_X{4,4}*full(Theta_obs_14'))...
                          + full(ghu(select_obs,:))*( full(XiXi_51)*full(Phi_obs_11')+full(XiXi_52)*full(Phi_obs_12')+full(XiXi_54)*full(Phi_obs_14')+full(XiXi_55)*full(Phi_obs_15'));
   % var-cov matrix of selected model variables, i.e., Gamma_0(y)                                                                                 
      Gamma_y_obs{1,1} = Gamma_obs{1,1} + Gamma_y_obs_second{1,1}+ Gamma_13_obs{1,1} + Gamma_13_obs{1,1}';                                                                               
                                                                                     
%--------------------------------------------------------------------------
% 4. Var-cov matrix between state vector and selected third order 
%    increments, i.e., Gamma_0(X,dy)
%--------------------------------------------------------------------------                                                                                     
temp = [ (full(Phi_X_11)*full(XiXi_11)+full(Phi_X_12)*full(XiXi_21)+full(Phi_X_14)*full(XiXi_41)+full(Phi_X_15)*full(XiXi_51))*full(Phi_obs_11')+(full(Phi_X_12)*full(XiXi_22)+full(Phi_X_14)*full(XiXi_42)+full(Phi_X_15)*full(XiXi_52)+full(Phi_X_11)*full(XiXi_21'))*full(Phi_obs_12')+(full(Phi_X_14)*full(XiXi_44)+full(Phi_X_15)*full(XiXi_54)+full(Phi_X_11)*full(XiXi_41')+full(Phi_X_12)*full(XiXi_42'))*full(Phi_obs_14') + (full(Phi_X_15)*full(XiXi_55) + full(Phi_X_11)*full(XiXi_51')+full(Phi_X_12)*full(XiXi_52')+full(Phi_X_14)*full(XiXi_54'))*full(Phi_obs_15')+full(Phi_X_13)*full(XiXi_33)*full(Phi_obs_13')
         (full(Phi_X_21)*full(XiXi_11)+full(Phi_X_22)*full(XiXi_21)+full(Phi_X_24)*full(XiXi_41))*full(Phi_obs_11')+(full(Phi_X_22)*full(XiXi_22)+full(Phi_X_24)*full(XiXi_42)+full(Phi_X_21)*full(XiXi_21'))*full(Phi_obs_12')+(full(Phi_X_24)*full(XiXi_44)+full(Phi_X_21)*full(XiXi_41')+full(Phi_X_22)*full(XiXi_42'))*full(Phi_obs_14')+(full(Phi_X_21)*full(XiXi_51')+full(Phi_X_22)*full(XiXi_52')+full(Phi_X_24)*full(XiXi_54'))*full(Phi_obs_15') + full(Phi_X_23)*full(XiXi_33)*full(Phi_obs_13')
         (full(Phi_X_31)*full(XiXi_11)+full(Phi_X_32)*full(XiXi_21))*full(Phi_obs_11')+(full(Phi_X_32)*full(XiXi_22)+full(Phi_X_31)*full(XiXi_21'))*full(Phi_obs_12')+(full(Phi_X_31)*full(XiXi_41')+full(Phi_X_32)*full(XiXi_42'))*full(Phi_obs_14')+(full(Phi_X_31)*full(XiXi_51')+full(Phi_X_32)*full(XiXi_52'))*full(Phi_obs_15')+ full(Phi_X_33)*full(XiXi_33)*full(Phi_obs_13')
         full(Phi_X_45)*full(XiXi_51)*full(Phi_obs_11')+full(Phi_X_45)*full(XiXi_52)*full(Phi_obs_12')+full(Phi_X_45)*full(XiXi_54)*full(Phi_obs_14')+full(Phi_X_45)*full(XiXi_55)*full(Phi_obs_15') ];

temp = [ (full(Theta_X_11)*Gamma_X{1,1}+full(Theta_X_12)*Gamma_X{2,1}+full(Theta_X_13)*Gamma_X{3,1}+full(Theta_X_14)*Gamma_X{4,1})*full(Theta_obs_11')+(full(Theta_X_12)*Gamma_X{2,2} + full(Theta_X_11)*Gamma_X{2,1}' + full(Theta_X_13)*Gamma_X{3,2}+full(Theta_X_14)*Gamma_X{4,2})*full(Theta_obs_12')+(full(Theta_X_13)*Gamma_X{3,3} + full(Theta_X_11)*Gamma_X{3,1}' + full(Theta_X_12)*Gamma_X{3,2}' + full(Theta_X_14)*Gamma_X{4,3})*full(Theta_obs_13')+(full(Theta_X_11)*Gamma_X{4,1}'+full(Theta_X_14)*Gamma_X{4,4} + full(Theta_X_12)*Gamma_X{4,2}' + full(Theta_X_13)*Gamma_X{4,3}')*full(Theta_obs_14')
         (full(Theta_X_22)*Gamma_X{2,1}+full(Theta_X_23)*Gamma_X{3,1}+full(Theta_X_24)*Gamma_X{4,1})*full(Theta_obs_11')+(full(Theta_X_22)*Gamma_X{2,2} + full(Theta_X_23)*Gamma_X{3,2} + full(Theta_X_24)*Gamma_X{4,2})*full(Theta_obs_12') + (full(Theta_X_23)*Gamma_X{3,3} + full(Theta_X_22)*Gamma_X{3,2}' + full(Theta_X_24)*Gamma_X{4,3})*full(Theta_obs_13') + (full(Theta_X_24)*Gamma_X{4,4} + full(Theta_X_22)*Gamma_X{4,2}' + full(Theta_X_23)*Gamma_X{4,3}')*full(Theta_obs_14')
         (full(Theta_X_33)*Gamma_X{3,1}+full(Theta_X_34)*Gamma_X{4,1})*full(Theta_obs_11')+(full(Theta_X_33)*Gamma_X{3,2}+full(Theta_X_34)*Gamma_X{4,2})*full(Theta_obs_12')+(full(Theta_X_33)*Gamma_X{3,3}+full(Theta_X_34)*Gamma_X{4,3})*full(Theta_obs_13')+(full(Theta_X_34)*Gamma_X{4,4} + full(Theta_X_33)*Gamma_X{4,3}')*full(Theta_obs_14')
         full(Theta_X_44)*(Gamma_X{4,1}*full(Theta_obs_11') + Gamma_X{4,2}*full(Theta_obs_12') + Gamma_X{4,3}*full(Theta_obs_13') + Gamma_X{4,4}*full(Theta_obs_14')) ]+temp;

Gamma_X_obs(:,1) = mat2cell(temp, [ns, ns^2, ns^3, ns], nobs); 

%--------------------------------------------------------------------------
% 5. Matrix of correlation and standard deviation of selected model
%    variables
%--------------------------------------------------------------------------
  % Standard deviation
    moments.third_order.standard_deviation = real( sqrt( diag( Gamma_y_obs{1,1} ) ) ) ;
  % Matrix of correlation
    moments.third_order.correlation = real( Gamma_y_obs{1,1}./(moments.third_order.standard_deviation*moments.third_order.standard_deviation') ); 
      
%--------------------------------------------------------------------------
% 6. Autocov matrices
%--------------------------------------------------------------------------     
    if options_.ar > 0     
        for j = 1: options_.ar
             Gamma_obs{1,j+1} = full(Theta_obs_11)*Gamma_X_obs{1,j}+full(Theta_obs_12)*Gamma_X_obs{2,j}+ full(Theta_obs_13)*Gamma_X_obs{3,j}+ full(Theta_obs_14)*Gamma_X_obs{4,j};
             Gamma_13_obs{1,j+1} = full(ghx(select_obs,:))* Gamma_X_obs{4,j} ;
             Gamma_X_obs{1,j+1} = full(Theta_X_11)*Gamma_X_obs{1,j}+full(Theta_X_12)*Gamma_X_obs{2,j}+full(Theta_X_13)*Gamma_X_obs{3,j}+full(Theta_X_14)*Gamma_X_obs{4,j};
             Gamma_X_obs{2,j+1} = full(Theta_X_22)*Gamma_X_obs{2,j}+full(Theta_X_23)*Gamma_X_obs{3,j}+full(Theta_X_24)*Gamma_X_obs{4,j};
             Gamma_X_obs{3,j+1} = full(Theta_X_33)*Gamma_X_obs{3,j}+full(Theta_X_34)*Gamma_X_obs{4,j};
             Gamma_X_obs{4,j+1} = full(Theta_X_44)*Gamma_X_obs{4,j};
             % autocov of selected model variables, i.e., Gamma_j(y)
               Gamma_y_obs{1,j+1} = Gamma_obs{1,j+1}+Gamma_y_obs_second{1,j+1}+Gamma_13_obs{1,j+1}+Gamma_13_obs{1,j+1}';
             % coefficients of autocorrelation of seleced model variables          
               moments.third_order.autocorrelation(:,j) = real( diag(  Gamma_y_obs{1,j+1}./(moments.third_order.standard_deviation*moments.third_order.standard_deviation')  ) );
        end
    end
    moments.third_order.Gamma_y_obs = Gamma_y_obs;
    
%--------------------------------------------------------------------------
% 7. Variance decomposition
%--------------------------------------------------------------------------
  
  % 7.1 Build coefficients
        % Theta_amp*A(D+)
          Theta_obs_amp_11 = ghx(select_obs,:)*(0.5*speye(ns));    
          Theta_obs_amp_12 = ghxx(select_obs,:);
          Theta_obs_amp_13 = (1/6)*ghxxx(select_obs,:);
          Theta_obs_amp_14 = 0.5*ghxuu(select_obs,:)*alt_kron(speye(ns),vec(Sigma_e));
        % Theta_risk*A(D+)  
          Theta_obs_risk_11 = ghx(select_obs,:)*(0.5*speye(ns));    
          Theta_obs_risk_14 = 0.5*ghxss_nlma(select_obs,:);
        % Phi_amp  
          Phi_obs_amp_11 = (1/6)*ghuuu(select_obs,:); 
          Phi_obs_amp_12 = 0.5*ghxxu(select_obs,:);
          Phi_obs_amp_13 = 0.5*ghxuu(select_obs,:);
          Phi_obs_amp_14 = ghxu(select_obs,:);        
  
  % 7.2. Amplification effect of third order increment, i.e., Gamma_0(amp)
         temp = (full(Phi_obs_amp_11)*full(XiXi_11)+full(Phi_obs_amp_12)*full(XiXi_21)+full(Phi_obs_amp_14)*full(XiXi_41))*full(Phi_obs_amp_11') + (full(Phi_obs_amp_12)*full(XiXi_22)+full(Phi_obs_amp_14)*full(XiXi_42)+full(Phi_obs_amp_11)*full(XiXi_21'))*full(Phi_obs_amp_12') + (full(Phi_obs_amp_14)*full(XiXi_44)+full(Phi_obs_amp_11)*full(XiXi_41')+full(Phi_obs_amp_12)*full(XiXi_42'))*full(Phi_obs_amp_14')+full(Phi_obs_amp_13)*full(XiXi_33)*full(Phi_obs_amp_13');
         Gamma_obs_amp{1,1} =  (full(Theta_obs_amp_11)*Gamma_X{1,1} + full(Theta_obs_amp_12)*Gamma_X{2,1} + full(Theta_obs_amp_13)*Gamma_X{3,1} + full(Theta_obs_amp_14)*Gamma_X{4,1})*full(Theta_obs_amp_11')...
                                +(full(Theta_obs_amp_12)*Gamma_X{2,2} + full(Theta_obs_amp_11)*Gamma_X{2,1}' + full(Theta_obs_amp_13)*Gamma_X{3,2} + full(Theta_obs_amp_14)*Gamma_X{4,2})*full(Theta_obs_amp_12')...
                                +(full(Theta_obs_amp_13)*Gamma_X{3,3} + full(Theta_obs_amp_11)*Gamma_X{3,1}' + full(Theta_obs_amp_12)*Gamma_X{3,2}' + full(Theta_obs_amp_14)*Gamma_X{4,3})*full(Theta_obs_amp_13')... 
                                +(full(Theta_obs_amp_11)*Gamma_X{4,1}' + full(Theta_obs_amp_14)*Gamma_X{4,4} + full(Theta_obs_amp_12)*Gamma_X{4,2}' + full(Theta_obs_amp_13)*Gamma_X{4,3}')*full(Theta_obs_amp_14')+temp; 
   
  % 7.3. Amplification effect of third order increment, i.e., Gamma_0(risk)
         temp = 0.5*full(ghuss_nlma(select_obs,:))*full(XiXi_55)*full(0.5*ghuss_nlma(select_obs,:)');
         Gamma_obs_risk{1,1} = (full(Theta_obs_risk_11)*Gamma_X{1,1} + full(Theta_obs_risk_14)*Gamma_X{4,1})*full(Theta_obs_risk_11') + (Theta_obs_risk_11*Gamma_X{4,1}' + Theta_obs_risk_14*Gamma_X{4,4})*full(Theta_obs_risk_14')+temp;

  % 7.4. Interplay between amplification and risk adjustment effect of
  %      third order increment, i.e., Gamma_0(amp,risk)
         Gamma_obs_amp_risk{1,1} = Gamma_obs{1,1} - Gamma_obs_amp{1,1} - Gamma_obs_risk{1,1};
  
  % 7.5. Amplification effect from interplay between first and third order 
  %      increments, i.e., Gamma_0((1)amp,(3)amp)
         Gamma_13_obs_amp{1,1} =  ghx(select_obs,:)*(   Gamma_X{4,1}*Theta_obs_amp_11' + Gamma_X{4,2}*Theta_obs_amp_12' + Gamma_X{4,3}*Theta_obs_amp_13' + Gamma_X{4,4}*Theta_obs_amp_14')...
                                  + ghu(select_obs,:)*( XiXi_51*Phi_obs_amp_11'+XiXi_52*Phi_obs_amp_12'+XiXi_54*Phi_obs_amp_14');
  
  % 7.6. Interplay between amplification and risk adjustment effect from 
  %      interaction between first and third order increments,
  %      i.e., Gamma_0((1)amp,(3)risk)
         Gamma_13_obs_amp_risk{1,1} = Gamma_13_obs{1,1} - Gamma_13_obs_amp{1,1};
  
  % 7.7. Total amplification effect
         Gamma_y_obs_amp{1,1} = Gamma_obs_amp{1,1}+ Gamma_y_obs_second{1,1} + Gamma_13_obs_amp{1,1} + Gamma_13_obs_amp{1,1}';
             
  % 7.8. Total interplay btw. amplification and risk adjustment effects
         Gamma_y_obs_amp_risk{1,1} = Gamma_obs_amp_risk{1,1}+ Gamma_13_obs_amp_risk{1,1}+ Gamma_13_obs_amp_risk{1,1}';
                                         
  % 7.9. Decomposition Results
         % Decomposition matrix
           %                                   Total risk adjustment effect  Total amplification effect    Total interplay btw. amp. and risk adjustment effect     
             moments.third_order.var_decomp = [ diag(Gamma_obs_risk{1,1})    diag(Gamma_y_obs_amp{1,1})    diag(Gamma_y_obs_amp_risk{1,1})        ];                       
           % Report results in percentage
             moments.third_order.var_decomp_report = (moments.third_order.var_decomp./repmat( full(diag(Gamma_y_obs{1,1})),[1,3])).*100;           
     
         % Complete decomposition matrix
           %                                             Total Risk Adj.            1st Order Amp.               2nd Order Incre. Amp.        3rd Order Incre. Amp.      1st-3rd Order Incre. Amp.                                3rd Order Incre. Interplay      1st-3rd Order Incre. Interplay
             moments.third_order.var_decomp_complete = [ diag(Gamma_obs_risk{1,1})  diag(Gamma_obs_first{1,1})   diag(Gamma_obs_second{1,1})  diag(Gamma_obs_amp{1,1})   diag(Gamma_13_obs_amp{1,1} + Gamma_13_obs_amp{1,1}')     diag(Gamma_obs_amp_risk{1,1})   diag(Gamma_13_obs_amp_risk{1,1}+Gamma_13_obs_amp_risk{1,1}') ];                
           % Report results in percentage
             moments.third_order.var_decomp_complete_report = (moments.third_order.var_decomp_complete./repmat( full(diag(Gamma_y_obs{1,1})),[1,7])).*100;

 

 warning on MATLAB:dividebyzero
        
        
        
        
        
        
        
        


   




