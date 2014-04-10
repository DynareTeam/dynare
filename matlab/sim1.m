function sim1()
% function sim1
% Performs deterministic simulations with lead or lag on one period.
% Uses sparse matrices.
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 1996-2013 Dynare Team
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

global M_ options_ oo_

endogenous_terminal_period = options_.endogenous_terminal_period;
vperiods = options_.periods*ones(1,options_.simul.maxit);
azero = options_.dynatol.f/1e7;

lead_lag_incidence = M_.lead_lag_incidence;

ny = M_.endo_nbr;

max_lag = M_.maximum_endo_lag;

nyp = nnz(lead_lag_incidence(1,:)) ;
iyp = find(lead_lag_incidence(1,:)>0) ;
ny0 = nnz(lead_lag_incidence(2,:)) ;
iy0 = find(lead_lag_incidence(2,:)>0) ;
nyf = nnz(lead_lag_incidence(3,:)) ;
iyf = find(lead_lag_incidence(3,:)>0) ;

nd = nyp+ny0+nyf;
nrc = nyf+1 ;
isp = [1:nyp] ;
is = [nyp+1:ny+nyp] ;
isf = iyf+nyp ;
isf1 = [nyp+ny+1:nyf+nyp+ny+1] ;
stop = 0 ;
iz = [1:ny+nyp+nyf];

periods = options_.periods;
steady_state = oo_.steady_state;
params = M_.params;
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;
i_cols_1 = nonzeros(lead_lag_incidence(2:3,:)');
i_cols_A1 = find(lead_lag_incidence(2:3,:)');
i_cols_T = nonzeros(lead_lag_incidence(1:2,:)');
i_cols_0 = nonzeros(lead_lag_incidence(2,:)');
i_cols_A0 = find(lead_lag_incidence(2,:)');
i_cols_j = 1:nd;
i_upd = ny+(1:periods*ny);

Y = endo_simul(:);

skipline()
disp (['-----------------------------------------------------']) ;
fprintf('MODEL SIMULATION:\n');

model_dynamic = str2func([M_.fname,'_dynamic']);
z = Y(find(lead_lag_incidence'));
[d1,jacobian] = model_dynamic(z,oo_.exo_simul, params, ...
                              steady_state,M_.maximum_lag+1);

A = sparse([],[],[],periods*ny,periods*ny,periods*nnz(jacobian));
res = zeros(periods*ny,1);

o_periods = periods;

ZERO = zeros(length(i_upd),1);

h1 = clock ;
for iter = 1:options_.simul.maxit
    h2 = clock ;
    
    i_rows = 1:ny;
    i_cols = find(lead_lag_incidence');
    i_cols_A = i_cols;

    for it = (M_.maximum_lag+1):(M_.maximum_lag+periods)

        [d1,jacobian] = model_dynamic(Y(i_cols), exo_simul, params, steady_state,it);
        if it == M_.maximum_lag+periods && it == M_.maximum_lag+1
            A(i_rows,i_cols_A0) = jacobian(:,i_cols_0);
        elseif it == M_.maximum_lag+periods
            A(i_rows,i_cols_A(i_cols_T)) = jacobian(:,i_cols_T);
        elseif it == M_.maximum_lag+1
            A(i_rows,i_cols_A1) = jacobian(:,i_cols_1);
        else
            A(i_rows,i_cols_A) = jacobian(:,i_cols_j);
        end

        res(i_rows) = d1;
        
        if endogenous_terminal_period && iter>1
            dr = max(abs(d1));
            if dr<azero
                vperiods(iter) = it;
                periods = it-M_.maximum_lag;
                break
            end
        end
    
        i_rows = i_rows + ny;
        i_cols = i_cols + ny;
        
        if it > M_.maximum_lag+1
            i_cols_A = i_cols_A + ny;
        end
    end
        
    err = max(abs(res));
    
    if options_.debug
        fprintf('\nLargest absolute residual at iteration %d: %10.3f\n',iter,err);    
        if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y))
            fprintf('\nWARNING: NaN or Inf detected in the residuals or endogenous variables.\n');    
        end
        skipline()
    end
    
    if err < options_.dynatol.f
        stop = 1 ;
        break
    end

    if endogenous_terminal_period && iter>1
        dy = ZERO;
        dy(1:i_rows(end)) = -A(1:i_rows(end),1:i_rows(end))\res(1:i_rows(end));
    else
        dy = -A\res;
    end
    
    Y(i_upd) =   Y(i_upd) + dy;

end

if endogenous_terminal_period
    err = evaluate_max_dynamic_residual(model_dynamic, Y, oo_.exo_simul, params, steady_state, o_periods, ny, max_lag, lead_lag_incidence);
    periods = o_periods;
end


if stop
    if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y))
        oo_.deterministic_simulation.status = 0;% NaN or Inf occurred
        oo_.deterministic_simulation.error = err;
        oo_.deterministic_simulation.iterations = iter;
        oo_.deterministic_simulation.periods = vperiods(1:iter);
        oo_.endo_simul = reshape(Y,ny,periods+2);
        skipline();
        fprintf('\nSimulation terminated after %d iterations.\n',iter);
        fprintf('Total time of simulation: %16.13f\n',etime(clock,h1));
        error('Simulation terminated with NaN or Inf in the residuals or endogenous variables. There is most likely something wrong with your model.');
    else
        skipline();
        fprintf('\nSimulation concluded successfully after %d iterations.\n',iter);
        fprintf('Total time of simulation: %16.13f\n',etime(clock,h1));
        fprintf('Max. Abs. Error         : %16.13f\n',err);
        fprintf('Convergency obtained!\n');
        oo_.deterministic_simulation.status = 1;% Convergency obtained.
        oo_.deterministic_simulation.error = err;
        oo_.deterministic_simulation.iterations = iter;
        oo_.deterministic_simulation.periods = vperiods(1:iter);
        oo_.endo_simul = reshape(Y,ny,periods+2);
    end
elseif ~stop
    skipline();
    fprintf('\nSimulation terminated after %d iterations.\n',iter);
    fprintf('Total time of simulation: %16.13f\n',etime(clock,h1));
    fprintf('Max. Abs. Error         : %16.13f\n',err);
    fprintf('WARNING : maximum number of iterations is reached (modify option maxit).\n') ;
    oo_.deterministic_simulation.status = 0;% more iterations are needed.
    oo_.deterministic_simulation.error = err;
    oo_.deterministic_simulation.periods = vperiods(1:iter);
    oo_.deterministic_simulation.iterations = options_.simul.maxit;
end
disp (['-----------------------------------------------------']) ;
skipline();