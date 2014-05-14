
function moments = nlma_th_moments(M_,oo_,options_,var_list_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nlma_th_momoments.m
%
% This is the main function that computes the theoretical moments, up to
% the third order accurate,  using nonlinear moving average solution.
% According to the order of approximation specified in options_.order, it
% calls
%  - nlma_th_mom_first.m
%  - nlma_th_mom_second.m
%  - nlma_th_mom_third.m
% to compute theoretical moments
%
% This file uses dyntable.m to report results
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

%--------------------------------------------------------------------------
% 1. Dynare version check and variable selectors
%--------------------------------------------------------------------------
  % Starting in Dynare 4.4.0, following fields are no longer in oo_.dr, 
  % they can be found in M_
    [ numeric_version ] = return_dynare_version( dynare_version );
    if numeric_version >= 4.4 
        oo_.dr.nstatic = M_.nstatic;
        oo_.dr.npred = M_.nspred; % note M_.nspred = M_.npred+M_.nboth;
        oo_.dr.nboth = M_.nboth;
        oo_.dr.nfwrd = M_.nfwrd;
    end
    
  % Build observables selector
    if size(var_list_,1) == 0
       var_list_ = M_.endo_names(1:M_.orig_endo_nbr, :);
    end
    nvar = size(var_list_,1);
    select_obs = zeros(nvar,1);
    for i=1:nvar
        i_tmp = strmatch(var_list_(i,:),M_.endo_names,'exact');
        if isempty(i_tmp)
           error (['One of the variables specified does not exist']) ;
        else
           select_obs(i) = i_tmp;
        end
    end
    ivar = select_obs; % For displaying nonstationary vars reminder
    select_obs = oo_.dr.inv_order_var(select_obs); 
    
  % Build state variable selector
    select_state = oo_.dr.nstatic + 1:oo_.dr.nstatic + oo_.dr.npred; 
    
%--------------------------------------------------------------------------
% 2. Save results for nlma moments calculation of all orders
%--------------------------------------------------------------------------
  % Save numbers of different type of variables
    moments.nstatic = oo_.dr.nstatic;
    moments.npred   = oo_.dr.npred;
    moments.nboth   = oo_.dr.nboth;
    moments.nfwrd   = oo_.dr.nfwrd;
    
  % Save variable selectors
    moments.select_state = select_state;
    moments.select_obs   = select_obs;
    
  % Number of state variables  
    moments.ns = size(select_state,2);
    
  % Number of observables
    moments.nobs = size(select_obs,1);
    
  % Number of shocks
    moments.ne = M_.exo_nbr;
  
%--------------------------------------------------------------------------
% 3. Compute nlma moments and index stationary variables
%--------------------------------------------------------------------------  
  if options_.order == 1
      moments = nlma_th_mom_first(moments,M_,oo_,options_);
      stationary_vars_first = find( diag(moments.first_order.Gamma_obs{1,1}) > 1e-12 );
      non_stationary_vars_first = setdiff(1:length(select_obs),stationary_vars_first);
  elseif options_.order == 2
      moments = nlma_th_mom_second(moments,M_,oo_,options_);
      stationary_vars_first = find( diag(moments.first_order.Gamma_obs{1,1}) > 1e-12 );
      stationary_vars_second = find( diag(moments.second_order.Gamma_y_obs{1,1}) > 1e-12 );
      non_stationary_vars_first = setdiff(1:length(select_obs),stationary_vars_first);
      non_stationary_vars_second = setdiff(1:length(select_obs),stationary_vars_second);
  elseif options_.order == 3
      if options_.pruning == 0;
          % This seperates the risk correction of the first order coefficients 
          % from the first order coefficients themselves and seperates out 
          % the blocks of the third order policy function ghxxx etc.
            nlma_oo = full_block_dr_new(oo_,M_,options_);
            moments.nlma_oo = nlma_oo; % Save for comparison with Dynare's oo_, could be removed
      else
          nlma_oo = oo_;
      end
      moments = nlma_th_mom_third(moments,M_,nlma_oo,options_);
      
      stationary_vars_first = find( diag(moments.first_order.Gamma_obs{1,1}) > 1e-12 );
      stationary_vars_second = find( diag(moments.second_order.Gamma_y_obs{1,1}) > 1e-12 );
      stationary_vars_third = find( diag(moments.third_order.Gamma_y_obs{1,1}) > 1e-12 );
      non_stationary_vars_first = setdiff(1:length(select_obs),stationary_vars_first);
      non_stationary_vars_second = setdiff(1:length(select_obs),stationary_vars_second);
      non_stationary_vars_third = setdiff(1:length(select_obs),stationary_vars_third);
  else
      error('This version only computes moments up to the third order of approximation.');
  end
  
  % Display nonstationary variable reminder
    if options_.order >=1 
               if length(non_stationary_vars_first) == length(select_obs)
                   disp( 'All variables appear either constant or nonstationary up to first order approximation ');
                   disp( 'Their first order accurate correlations and auto-correlations will not be displayed');
               elseif isempty(non_stationary_vars_first) == 0 && length(non_stationary_vars_first) ~= length(select_obs)
                   disp( 'Following variables appear either constant or nonstationary up to first order approximation ');
                   disp( M_.endo_names(ivar(non_stationary_vars_first,:),:) );
                   disp( 'Their first order accurate correlations and auto-correlations will not be displayed');      
               end
    end
    if options_.order >=2
               if length(non_stationary_vars_second) == length(select_obs)
                   disp( 'All variables appear either constant or nonstationary up to second order approximation ');
                   disp( 'Their second order accurate correlations and auto-correlations will not be displayed');
               elseif isempty(non_stationary_vars_second) == 0 && length(non_stationary_vars_second) ~= length(select_obs)
                   disp( 'Following variables appear either constant or nonstationary up to second order approximation ');
                   disp( M_.endo_names(ivar(non_stationary_vars_second,:),:) );
                   disp( 'Their second order accurate correlations and auto-correlations will not be displayed');      
               end
    end
    if options_.order == 3
               if length(non_stationary_vars_third) == length(select_obs)
                   disp( 'All variables appear either constant or nonstationary up to third order approximation ');
                   disp( 'Their third order accurate correlations and auto-correlations will not be displayed');
               elseif isempty(non_stationary_vars_third) == 0 && length(non_stationary_vars_third) ~= length(select_obs)
                   disp( 'Following variables appear either constant or nonstationary up to third order approximation ');
                   disp( M_.endo_names(ivar(non_stationary_vars_third,:),:) );
                   disp( 'Their third order accurate correlations and auto-correlations will not be displayed');      
               end
    end
    
%--------------------------------------------------------------------------
% 4. Display results
%--------------------------------------------------------------------------  
  % Table : mean
    if options_.order == 1
         disp_mean = moments.first_order.mean;
         headers = char('VARIABLE','MEAN(1st Order Accurate)');
    elseif options_.order == 2 || options_.order == 3
         disp_mean = [ moments.first_order.mean moments.second_order.mean ];
         headers = char('VARIABLE','MEAN(1st Order Accurate)','MEAN(2nd and 3rd Order Accurate)');
    end
    title   = 'NLMA THEORETICAL MEAN'; 
    labels  = deblank( M_.endo_names(oo_.dr.order_var(select_obs),:) );
    lh      = size(labels,2)+2;
    dyntable( title, headers, labels, disp_mean , lh, 11, 4 );
    
  % Table : mean decomposition
    if options_.order == 2 || options_.order == 3
        title   = 'NLMA MEAN(2nd and 3rd Order Accurate) DECOMPOSITION IN LEVEL';
        headers = char('VARIABLE','Mean(2nd and 3rd Order Accurate)', ...
                                  'Det. Steady State',...
                                  'Risk Adj. from Var. of Future Shock',...
                                  'Risk Adj. from Var. of Past Shock'       );
        labels  = deblank( M_.endo_names(oo_.dr.order_var(select_obs),:) );
        lh      = size(labels,2)+2;
        dyntable( title, headers, labels, moments.second_order.mean_decomp, lh, 11, 4 );
    end
 
  % Table : standard deviation and variance
    if options_.order == 1
        disp_std_var = [ moments.first_order.standard_deviation  diag(moments.first_order.Gamma_obs{1,1}) ];
        headers = char('VARIABLE','STD. (1st Order Accurate.)','VAR.(1st Order Accurate)');
    elseif options_.order == 2
        disp_std_var = [ moments.first_order.standard_deviation    moments.second_order.standard_deviation ...
                         diag(moments.first_order.Gamma_obs{1,1})  diag(moments.second_order.Gamma_y_obs{1,1})    ];
        headers = char( 'VARIABLE' , 'STD.(1st Order Accurate)' , 'STD.(2nd Order Accurate)' , ...
                                     'VAR.(1st Order Accurate)' , 'VAR.(2nd Order Accurate)'      );
    elseif options_.order == 3
        disp_std_var = [ moments.first_order.standard_deviation    moments.second_order.standard_deviation     moments.third_order.standard_deviation ...
                         diag(moments.first_order.Gamma_obs{1,1})  diag(moments.second_order.Gamma_y_obs{1,1}) diag(moments.third_order.Gamma_y_obs{1,1})   ];
                            
        headers = char( 'VARIABLE' , 'STD.(1st Order Accurate)' , 'STD.(2nd Order Accurate)' , 'STD.(3rd Order Accurate)' , ...
                                     'VAR.(1st Order Accurate)' , 'VAR.(2nd Order Accurate)' , 'VAR.(3rd Order Accurate)'      );
                                     
    end
    title   = 'NLMA THEORETICAL STANDARD DEVIATION AND VARIANCE';
    labels  = deblank( M_.endo_names(oo_.dr.order_var(select_obs),:) );
    lh      = size(labels,2)+2;
    dyntable( title, headers, labels, disp_std_var , lh, 11, 4 );
    
  % Table : correlation
    if size( stationary_vars_first,1 ) > 0 && options_.order >= 1                
        title   = 'NLMA THEORETICAL CORRELATIONS(1st Order Accurate)';
        labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_first)),:));
        headers = char('Variables',labels);
        dyntable( title, headers, labels, ...
                  moments.first_order.correlation(stationary_vars_first,stationary_vars_first), ...
                  lh, 8, 4 );
    end
    if options_.order >= 2
        if size( stationary_vars_second,1 ) > 0 
           title   = 'NLMA THEORETICAL CORRELATIONS(2nd Order Accurate)';
           labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_second)),:));
           headers = char('Variables',labels);
           dyntable( title, headers, labels, ...
                     moments.second_order.correlation(stationary_vars_second,stationary_vars_second), ...
                     lh, 8, 4 );
        end
    end
    if options_.order == 3
        if size( stationary_vars_third,1 ) > 0
           title   = 'NLMA THEORETICAL CORRELATIONS(3rd Order Accurate)';
           labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_third)),:));
           headers = char('Variables',labels);
           dyntable( title, headers, labels, ...
                     moments.third_order.correlation(stationary_vars_third,stationary_vars_third), ...
                     lh, 8, 4 );
        end
    end
 
  % Table : coefficient of autocorrelation
    if options_.ar > 0          
         if size( stationary_vars_first,1 ) > 0 && options_.order >= 1
             title   = 'NLMA THEORETICAL COEFFICIENTS OF AUTOCORRELATION(1st Order Accurate) ';
             headers = char( 'Order ',int2str([1:options_.ar]'));
             labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_first)),:));
             lh      = size(labels,2)+2;
             dyntable( title, headers, labels, moments.first_order.autocorrelation(stationary_vars_first,:), lh, 8, 4 );
         end
         if options_.order >= 2
             if size( stationary_vars_second,1 ) > 0 
                title   = 'NLMA THEORETICAL COEFFICIENTS OF AUTOCORRELATION(2nd Order Accurate) ';
                headers = char( 'Order ',int2str([1:options_.ar]'));
                labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_second)),:));
                lh      = size(labels,2)+2;
                dyntable( title, headers, labels, moments.second_order.autocorrelation(stationary_vars_second,:), lh, 8, 4 );
             end
         end
         if options_.order == 3
             if size( stationary_vars_third,1 ) > 0 
                 title   = 'NLMA THEORETICAL COEFFICIENTS OF AUTOCORRELATION(3rd Order Accurate) ';
                 headers = char( 'Order ',int2str([1:options_.ar]'));
                 labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_third)),:));
                 lh      = size(labels,2)+2;
                 dyntable( title, headers, labels, moments.third_order.autocorrelation(stationary_vars_third,:), lh, 8, 4 );
             end
         end
    end
    
  % Table : variance decomposition
    if  options_.order == 3
        if size( stationary_vars_third,1 ) > 0 
        % Decomposition 
          title   = 'NLMA VARIANCE (3rd Order Accurate) DECOMPOSITION IN PERCENTAGE';
          labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_third)),:));
          headers = char( 'VARIABLE' , 'Total Risk Adj.', 'Total Amp.', 'Total Risk-Amp. Interplay');
          dyntable( title, headers, labels, moments.third_order.var_decomp_report(stationary_vars_third,:) ,lh, 8, 2 );
        % Complete decomposition 
          title   = 'NLMA COMPLETE VARIANCE (3rd Order Accurate) DECOMPOSITION IN PERCENTAGE';
          labels  = deblank(M_.endo_names(oo_.dr.order_var(select_obs(stationary_vars_third)),:));
          headers = char( 'VARIABLE' , 'Total Risk Adj.','1st Order Amp.','2nd Order Incre. Amp.','3rd Order Incre. Amp.','1st-3rd Order Incre. Amp.','3rd Order Incre. Risk-Amp. Interplay','1st-3rd Order Incre. Risk-Amp. Interplay');
          dyntable( title, headers, labels, moments.third_order.var_decomp_complete_report(stationary_vars_third,:) ,lh, 8, 2 );
        end
    end








  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  