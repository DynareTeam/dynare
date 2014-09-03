
function nlma_simul = nlma_simul(M_,options_,var_list_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nlma_simul_plot.m
%
% This file computes and plots the nonlinear moving average simulation, 
% including a decomposition of the simulation into contributing orders 
% of nonlinearity and uncertainty.
%
% According to the order of approximation specified in options_.order, it
% calls
%  - pruning_abounds.m
% to compute nlma simulation
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
% 2. Compute nlma simulations
%--------------------------------------------------------------------------
  % Length of simulation
    T    = options_.periods;
    drop = options_.drop;
   
  % Construct x-axis for all simulation plot
    simul_xaxis = 1+drop : T; 
  
  % Import shock sequence for simulation from Dynare
    simul_shock_sequence = oo_.exo_simul';
    
  % Compute nlma simulation
    simul_data = pruning_abounds( M_, options_, simul_shock_sequence, T, options_.order, 'lan_meyer-gohde' );
    nlma_simul.simul_data = simul_data;

%--------------------------------------------------------------------------
% 3. Plot results
%--------------------------------------------------------------------------

 if options_.nograph == 0
     
    for ii = 1:length(variable_select) % Loop over all selected variables
          
          % Plot first order accuate simulation
            if options_.order == 1
               figure('Units','characters','Position',[1 1 240 70]);
                  plot( simul_xaxis, simul_data.total(variable_select(ii),drop+1:end), 'b' );
                  eval(sprintf('title(''Simulation of %s'')',M_.endo_names(variable_select(ii),:)))
                  ylabel('Value');  xlabel('Periods'); 
                  %axis tight;
                  xlim([ 1+drop T]);
                  legend({'First-Order Accurate'},'Location','Best')

              
          % Plot first and second order accurate simulation and its
          % decomposition, and plot Dynare's (non)pruned simulation for
          % comparison
            elseif options_.order==2
              figure('Units','characters','Position',[1 1 240 70]);
                % First and second order accurate simulations
                  subplot(3,3,[1:6]); 
                       plot( simul_xaxis, simul_data.total( variable_select(ii),drop+1:end ),'b',...      % Second order accurate simulation
                             simul_xaxis, simul_data.first( variable_select(ii),drop+1:end )     ...      % First order accurate simulation
                                          + oo_.dr.ys(variable_select(ii),:),'b-.',              ...
                             simul_xaxis, oo_.endo_simul(variable_select(ii),drop+1:end),'k--'        );  % Dynare's (non)pruned second order accurate simulation
                       eval(sprintf('title(''Simulation of %s'')',M_.endo_names(variable_select(ii),:)));
                       ylabel('Value');  xlabel('Periods'); 
                       %axis tight;
                       xlim([ 1+drop T]);
                       if options_.pruning == 1
                           legend({'Second-Order Accurate','First-Order Accurate','Dynare: Pruned Second Order'},'Location','Best');
                       else
                           legend({'Second-Order Accurate','First-Order Accurate','Dynare: Second Order'},'Location','Best');
                       end
                % First order component of second order accurate simulation
                  subplot(3,3,[7]); 
                       plot( simul_xaxis, simul_data.first(variable_select(ii),drop+1:end),'b')
                       title('First-Order Component'); ylabel('Deviations'); xlabel('Periods');
                       hold on;
                       plot( simul_xaxis,zeros(1, T-drop),'k-');
                       hold off; 
                       %axis tight;
                       xlim([ 1+drop T]);
                % Second order component of second order accurate simulation
                  subplot(3,3,[8]); 
                       plot( simul_xaxis, simul_data.second(variable_select(ii),drop+1:end),'b')
                       title('Second-Order Component'); ylabel('Deviations');  xlabel('Periods');
                       hold on;
                       plot( simul_xaxis, zeros(1, T-drop),'k-');
                       hold off; 
                       %axis tight;
                       xlim([ 1+drop T]);
                % Constant risk correction term
                  ghs2_nlma = oo_.dr.ghs2_nlma(oo_.dr.inv_order_var,:);
                  subplot(3,3,[9]); 
                       plot( simul_xaxis, repmat(oo_.dr.ys(variable_select(ii)), [1 T-drop]),'b',...
                             simul_xaxis, repmat(oo_.dr.ys(variable_select(ii)) + 0.5*ghs2_nlma(variable_select(ii)), [1 T-drop]),'b-.');
                       title('Steady-State and Constant'); 
                       legend({'Steady-State','plus Risk Adjutsment'});
                       xlim([ 1+drop T]);
                
          % Plot first, second and third order accurate simulations and its
          % decomposition, and plot Dynare's (non)pruned simulation for 
          % comparison.
            elseif options_.order==3
              figure('Units','characters','Position',[1 1 240 70]);
                % First, second and third order accurate simulations
                  ghs2_nlma = oo_.dr.ghs2_nlma(oo_.dr.inv_order_var,:);
                  subplot( 4,3,[1:6]); 
                       plot( simul_xaxis, simul_data.total(variable_select(ii),drop+1:end),'b', ...          % Third order accurate simulation
                             simul_xaxis, simul_data.first(variable_select(ii),drop+1:end) ...               % Second order accurate simulation
                                          + simul_data.second(variable_select(ii),drop+1:end)...
                                          + oo_.dr.ys(variable_select(ii),:)...
                                          + 0.5*ghs2_nlma(variable_select(ii),:),'b--',...
                             simul_xaxis, simul_data.first(variable_select(ii),drop+1:end)...                % First order accurate simulation
                                          + oo_.dr.ys(variable_select(ii),:),'b-.',...
                             simul_xaxis, oo_.endo_simul(variable_select(ii),drop+1:end),'k--'           );  % Dynare's (non)pruned third order accurate simulation
                       eval(sprintf('title(''Simulation of %s'')',M_.endo_names(variable_select(ii),:)));
                       ylabel('Value');  xlabel('Periods');
                       if options_.pruning == 1
                           legend({'Third-Order Accurate','Second-Order Accurate','First-Order Accurate','Dynare: Pruned Third Order'},'Location','Best');
                       else
                           legend({'Third-Order Accurate','Second-Order Accurate','First-Order Accurate','Dynare: Third Order'},'Location','Best');
                       end
                       hold on;
                       plot( simul_xaxis, zeros(1, T-drop),'k-');
                       hold off; 
                       %axis tight;
                       xlim([ 1+drop T]);
                % First order component of third order accurate simulation
                  subplot(4,3,[7]); 
                       plot( simul_xaxis, simul_data.first(variable_select(ii),drop+1:end),'b');
                       title('First-Order Component');  ylabel('Deviations');   xlabel('Periods');
                       hold on;
                       plot( simul_xaxis, zeros(1, T-drop),'k-');
                       hold off; 
                       %axis tight;
                       xlim([ 1+drop T]);
                % Second order component of third order accurate simulation
                  subplot(4,3,[8]); 
                       plot( simul_xaxis, simul_data.second(variable_select(ii),drop+1:end),'b')
                       title('Second-Order Component');    ylabel('Deviations');   xlabel('Periods');
                       hold on;
                       plot( simul_xaxis, zeros(1, T-drop),'k-');
                       hold off; 
                       %axis tight;
                       xlim([ 1+drop T]);
                % Time varying risk correct term
                  subplot(4,3,[10]); 
                       plot( simul_xaxis,simul_data.first_sigma_2(variable_select(ii),drop+1:end),'b');
                       title('Risk Correction to First-Order'); ylabel('Deviations'); xlabel('Periods');
                       hold on;
                       plot( simul_xaxis,zeros(1, T-drop),'k-');
                       hold off; 
                       %axis tight;
                       xlim([ 1+drop T]);
                % Third order component of third order accurate simulation
                  subplot(4,3,[11]); 
                       plot( simul_xaxis, simul_data.third(variable_select(ii),drop+1:end),'b');
                       title('Third-Order Component');  ylabel('Deviations');   xlabel('Periods');
                       hold on;
                       plot( simul_xaxis, zeros(1, T-drop),'k-'); 
                       %axis tight;
                       xlim([ 1+drop T]);
                       hold off;
                % Constant risk correction term
                  ghs2_nlma = oo_.dr.ghs2_nlma(oo_.dr.inv_order_var,:);
                  subplot(4,3,[9]); 
                       plot( simul_xaxis, repmat(oo_.dr.ys(variable_select(ii)), [1 T-drop]),'b',...
                             simul_xaxis, repmat(oo_.dr.ys(variable_select(ii))+ 0.5*ghs2_nlma(variable_select(ii)), [1 T-drop]),'b-.');
                       title('Steady-State and Constant');
                       legend({'Steady-State','plus Risk Adjustment'});
                       xlim([ 1+drop T]);
            end
          
     end  % Loop over variables ends
    
 end % nograph check ends
 
  
