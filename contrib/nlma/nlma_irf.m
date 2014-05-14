
function nlma_irf = nlma_irf(M_,options_,var_list_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nlma_irf_plot.m
%
% This file computes and plots the nonlinear moving average impulse 
% response, including a decomposition of the IRF into contributing orders 
% of nonlinearity and uncertainty.
%
% According to the order of approximation specified in options_.order, it
% calls
%  - pruning_abounds.m
% to compute nlma irfs
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
global oo_
%--------------------------------------------------------------------------
% 1. Build variable selector
%--------------------------------------------------------------------------
  if isempty(var_list_)==0
     [waste, variable_select] = ismember(cellstr(var_list_), cellstr(M_.endo_names));
  else
     variable_select=1:M_.endo_nbr;
  end
    
%--------------------------------------------------------------------------
% 2. Compute nlma impulse response
%--------------------------------------------------------------------------
  % Length of irf
    T = options_.irf;
   
  % Scale size of the shock    
    irf_shock_scale = 1;
    
  % Construct x-axis for all irf plot
    irf_xaxis = 1 : T; 
  
  % Preallocation
    irf_all = cell(M_.exo_nbr,1);
    
%   % Compute irf
%     for jj = 1 : M_.exo_nbr
%         irf_shock_sequence = zeros(M_.exo_nbr,T); % Pre-allocate and reset irf shock sequence
%         irf_shock_sequence(jj,1) = irf_shock_scale*(M_.Sigma_e(jj,jj))^(1/2);
%         irf_all{jj,1} = pruning_abounds( M_, options_, irf_shock_sequence, T, options_.order, 'lan_meyer-gohde' );
%     end
%     nlma_irf.irf_all = irf_all;
    
  % Compute irf, allowing correlated shocks
    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord) = M_.Sigma_e+1e-14*eye(M_.exo_nbr);
    cs = transpose(chol(SS));
    irf_shocks_indx = getIrfShocksIndx();
    for jj = irf_shocks_indx
        irf_shock_sequence = zeros(M_.exo_nbr,T); % Pre-allocate and reset irf shock sequence
        irf_shock_sequence(:,1) = irf_shock_scale*cs(M_.exo_names_orig_ord,jj);
        irf_all{jj,1} = pruning_abounds( M_, options_, irf_shock_sequence, T, options_.order, 'lan_meyer-gohde' );
    end
    nlma_irf.irf_all = irf_all;    
    
%--------------------------------------------------------------------------
% 3. Plot results
%--------------------------------------------------------------------------
  if options_.nograph == 0

    % for jj = 1:M_.exo_nbr % Loop over shocks
    for jj = irf_shocks_indx % Loop over (maybe correlated) shocks
      
      irf = irf_all{jj,1};
      
      for ii = 1:length(variable_select) % For each shock, loop over all selected variables
          
          % Plot first order accuate irf
            if options_.order == 1
               figure('Units','characters','Position',[1 1 240 70]);
                 plot( irf_xaxis, irf.first(variable_select(ii),:), 'b' );
                 eval( sprintf('title(''Impulse Response of %s to a %s Std. Dev. Shock in %s'')',...
                       M_.endo_names(variable_select(ii),:),num2str((irf_shock_scale)),M_.exo_names(jj,:)) );
                 ylabel('Deviations');  xlabel('Periods since Shock Realization');
                 legend({'First-Order Accurate'});
                 hold on;
                 plot(irf_xaxis,zeros(1, T),'k-');
                 hold off; 
                 %axis tight;
                 xlim([1 T]);
              
          % Plot second order accurate irf and its decomposition
            elseif options_.order==2
              figure('Units','characters','Position',[1 1 240 70]);
                % First and second order accurate irfs
                  subplot(3,3,[1:6]); 
                       plot( irf_xaxis, irf.first(variable_select(ii),:)+irf.second(variable_select(ii),:),'b',...
                             irf_xaxis, irf.first(variable_select(ii),:),'b-.');
                       eval( sprintf('title(''Impulse Response of %s to a %s Std. Dev. Shock in %s'')',...
                             M_.endo_names(variable_select(ii),:),num2str((irf_shock_scale)),M_.exo_names(jj,:)));
                       ylabel('Deviations');  xlabel('Periods since Shock Realization');
                       legend({'Second-Order Accurate','First-Order Accurate'})
                       hold on;
                       plot( irf_xaxis,zeros(1, T),'k-' );
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % First order component of second order accurate irf
                  subplot(3,3,[7]); 
                       plot(irf_xaxis,irf.first(variable_select(ii),:),'b')
                       title('First-Order Component'); ylabel('Deviations'); xlabel('Periods');
                       hold on;
                       plot(irf_xaxis,zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % Second order component of second order accurate irf
                  subplot(3,3,[8]); 
                       plot(irf_xaxis,irf.second(variable_select(ii),:),'b')
                       title('Second-Order Component'); ylabel('Deviations');  xlabel('Periods');
                       hold on;
                       plot(irf_xaxis,zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % Constant risk correction term
                  ghs2_nlma = oo_.dr.ghs2_nlma(oo_.dr.inv_order_var,:);
                  subplot(3,3,[9]); 
                       plot( irf_xaxis, repmat(oo_.dr.ys(variable_select(ii)), [1 T]),'b',...
                             irf_xaxis, repmat(oo_.dr.ys(variable_select(ii)) + 0.5*ghs2_nlma(variable_select(ii)), [1 T]),'b-.');
                       title('Steady-State and Constant'); 
                       legend({'Steady-State','plus Risk Adjutsment'});
                       xlim([1 T]);
                
          % Plot third order accurate irf and its decomposition
            elseif options_.order==3
              figure('Units','characters','Position',[1 1 240 70]);
                % First, second and third order accurate irfs
                  subplot( 4,3,[1:6]); 
                       plot( irf_xaxis, irf.first(variable_select(ii),:)   ...
                                        + irf.second(variable_select(ii),:)...
                                        + irf.third(variable_select(ii),:) ...
                                        + irf.first_sigma_2(variable_select(ii),:),'b', ...
                             irf_xaxis, irf.first(variable_select(ii),:)   ...
                                        + irf.second(variable_select(ii),:),'b--',...
                             irf_xaxis, irf.first(variable_select(ii),:),'b-.'         );
                       eval( sprintf('title(''Impulse Response of %s to a %s Std. Dev. Shock in %s'')',...
                       M_.endo_names(variable_select(ii),:),num2str((irf_shock_scale)),M_.exo_names(jj,:)) );
                       ylabel('Deviations');   xlabel('Periods since Shock Realization');
                       legend({'Third-Order Accurate','Second-Order Accurate','First-Order Accurate'})
                       hold on;
                       plot( irf_xaxis,zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % First order component of third order accurate irf
                  subplot(4,3,[7]); 
                       plot( irf_xaxis, irf.first(variable_select(ii),:),'b');
                       title('First-Order Component');  ylabel('Deviations');   xlabel('Periods');
                       hold on;
                       plot( irf_xaxis, zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % Second order component of third order accurate irf
                  subplot(4,3,[8]); 
                       plot( irf_xaxis, irf.second(variable_select(ii),:),'b')
                       title('Second-Order Component');    ylabel('Deviations');   xlabel('Periods');
                       hold on;
                       plot( irf_xaxis, zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % Time varying risk correct term
                  subplot(4,3,[10]); 
                       plot( irf_xaxis,irf.first_sigma_2(variable_select(ii),:),'b');
                       title('Risk Correction to First-Order'); ylabel('Deviations'); xlabel('Periods');
                       hold on;
                       plot( irf_xaxis,zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % Third order component of third order accurate irf
                  subplot(4,3,[11]); 
                       plot( irf_xaxis, irf.third(variable_select(ii),:),'b');
                       title('Third-Order Component');  ylabel('Deviations');   xlabel('Periods');
                       hold on;
                       plot( irf_xaxis, zeros(1, T),'k-');
                       hold off; 
                       %axis tight;
                       xlim([1 T]);
                % Constant risk correction term
                  ghs2_nlma = oo_.dr.ghs2_nlma(oo_.dr.inv_order_var,:);
                  subplot(4,3,[9]); 
                       plot( irf_xaxis, repmat(oo_.dr.ys(variable_select(ii)), [1 T]),'b',...
                             irf_xaxis, repmat(oo_.dr.ys(variable_select(ii))+ 0.5*ghs2_nlma(variable_select(ii)), [1 T]),'b-.');
                       title('Steady-State and Constant');
                       legend({'Steady-State','plus Risk Adjustment'}); 
                       xlim([1 T]);
            end
          
      end  % Loop over variables ends
    end % Loop over shocks ends
  
  end