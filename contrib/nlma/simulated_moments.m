
function simul_moments = simulated_moments(M_,oo_, options_,var_list_,pruning_order,pruning_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulated_moments.m
%
% This file simulates the chosen pruning algorithm (specified by the
% option, pruning_type, see the documentation of pruning_abounds.m for
% available pruning types), and computes moments using the simulated data
%
% According to the order of approximation specified in options_.order, it
% calls
%  - pruning_abounds.m
% to compute simulation
%
% The moment calculation part of this file is taken from disp_moments.m
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
% Copyright (C) 2001-2012 Dynare Team
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
% 1. Build variable selector
%--------------------------------------------------------------------------
  if size(var_list_,1) == 0
      var_list_ = M_.endo_names(1:M_.orig_endo_nbr, :);
  end
  nvar = size(var_list_,1);
  ivar = zeros(nvar,1);
  for i=1:nvar
      i_tmp = strmatch(var_list_(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
          error (['One of the variable specified does not exist']) ;
      else
          ivar(i) = i_tmp;
      end
  end

%--------------------------------------------------------------------------
% 2. Simulation
%--------------------------------------------------------------------------
  shock_sequence = oo_.exo_simul';
  simul_length   = options_.periods;
  y = pruning_abounds(M_,options_,shock_sequence, simul_length,pruning_order,pruning_type);
  y = y.total(ivar,options_.drop+1:end)';

%--------------------------------------------------------------------------
% 3. Compute moments
%--------------------------------------------------------------------------
  m = mean(y);

  if options_.hp_filter
     [hptrend,y] = sample_hp_filter(y,options_.hp_filter);
  else
      y = bsxfun(@minus, y, m);
  end

  s2 = mean(y.*y);
  s = sqrt(s2);
  simul_moments.mean = transpose(m);
  simul_moments.var = y'*y/size(y,1);

  labels = deblank(M_.endo_names(ivar,:));
  if options_.nomoments == 0
      z = [ m' s' s2' (mean(y.^3)./s2.^1.5)' (mean(y.^4)./(s2.*s2)-3)' ];    
      
      title = sprintf('MOMENTS OF SIMULATED VARIABLES, %s, ORDER %s', pruning_type, num2str(pruning_order) ) ;
      if options_.hp_filter
          title = [title ' (HP filter, lambda = ' ...
                   num2str(options_.hp_filter) ')'];
      end
      headers=char('VARIABLE','MEAN','STD. DEV.','VARIANCE','SKEWNESS', ...
                   'KURTOSIS');
      dyntable(title,headers,labels,z,size(labels,2)+2,16,6);
  end

  if options_.nocorr == 0
     corr = (y'*y/size(y,1))./(s'*s);
     if options_.noprint == 0
        title = sprintf('CORRELATION OF SIMULATED VARIABLES, %s, ORDER %s', pruning_type, num2str(pruning_order)  ) ;
        if options_.hp_filter
            title = [title ' (HP filter, lambda = ' ...
                     num2str(options_.hp_filter) ')'];
        end
        headers = char('VARIABLE',M_.endo_names(ivar,:));
        dyntable(title,headers,labels,corr,size(labels,2)+2,8,4);
     end
  end

  ar = options_.ar;
  if ar > 0
     autocorr = [];
     for i=1:ar
        simul_moments.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));
        autocorr = [ autocorr diag(simul_moments.autocorr{i}) ];
     end
     if options_.noprint == 0
        title = sprintf('AUTOCORRELATION OF SIMULATED VARIABLES, %s, ORDER %s', pruning_type, num2str(pruning_order)  ) ;
        if options_.hp_filter
            title = [title ' (HP filter, lambda = ' ...
                     num2str(options_.hp_filter) ')'];
        end
        headers = char('VARIABLE',int2str([1:ar]'));
        dyntable(title,headers,labels,autocorr,size(labels,2)+2,8,4);
     end
  end

