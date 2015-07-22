function oo_ = sim1(options_, M_, oo_)
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

% Copyright (C) 1996-2015 Dynare Team
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

verbose = options_.verbosity;

endogenous_terminal_period = options_.endogenous_terminal_period;
vperiods = options_.periods*ones(1,options_.simul.maxit);
azero = options_.dynatol.f/1e7;

lead_lag_incidence = M_.lead_lag_incidence;

ny = M_.endo_nbr;

maximum_lag = M_.maximum_lag;
max_lag = M_.maximum_endo_lag;

nyp = nnz(lead_lag_incidence(1,:)) ;
ny0 = nnz(lead_lag_incidence(2,:)) ;
nyf = nnz(lead_lag_incidence(3,:)) ;

nd = nyp+ny0+nyf;
stop = 0 ;

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
i_cols_j = (1:nd)';
i_upd = maximum_lag*ny+(1:periods*ny);

Y = endo_simul(:);

if verbose
    skipline()
    printline(56)
    disp('MODEL SIMULATION:')
    skipline()
end

model_dynamic = str2func([M_.fname,'_dynamic']);
z = Y(find(lead_lag_incidence'));
[d1,jacobian] = model_dynamic(z,oo_.exo_simul, params, ...
                              steady_state,maximum_lag+1);

res = zeros(periods*ny,1);

o_periods = periods;

ZERO = zeros(length(i_upd),1);

h1 = clock ;
iA = zeros(periods*M_.NNZDerivatives(1),3);
for iter = 1:options_.simul.maxit
    h2 = clock ;

    i_rows = (1:ny)';
    i_cols_A = find(lead_lag_incidence');
    i_cols = i_cols_A+(maximum_lag-1)*ny;
    m = 0;
    for it = (maximum_lag+1):(maximum_lag+periods)
        [d1,jacobian] = model_dynamic(Y(i_cols), exo_simul, params, ...
                                      steady_state,it);
        if it == maximum_lag+periods && it == maximum_lag+1
            [r,c,v] = find(jacobian(:,i_cols_0));
            iA((1:length(v))+m,:) = [i_rows(r),i_cols_A0(c),v];
        elseif it == maximum_lag+periods
            [r,c,v] = find(jacobian(:,i_cols_T));
            iA((1:length(v))+m,:) = [i_rows(r),i_cols_A(i_cols_T(c)),v];
        elseif it == maximum_lag+1
            [r,c,v] = find(jacobian(:,i_cols_1));
            iA((1:length(v))+m,:) = [i_rows(r),i_cols_A1(c),v];
        else
            [r,c,v] = find(jacobian(:,i_cols_j));
            iA((1:length(v))+m,:) = [i_rows(r),i_cols_A(c),v];
        end
        m = m + length(v);

        res(i_rows) = d1;

        if endogenous_terminal_period && iter>1
            dr = max(abs(d1));
            if dr<azero
                vperiods(iter) = it;
                periods = it-maximum_lag+1;
                break
            end
        end

        i_rows = i_rows + ny;
        i_cols = i_cols + ny;

        if it > maximum_lag+1
            i_cols_A = i_cols_A + ny;
        end
    end

    err = max(abs(res));

    if options_.debug
        fprintf('\nLargest absolute residual at iteration %d: %10.3f\n',iter,err);
        if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y))
            fprintf('\nWARNING: NaN or Inf detected in the residuals or endogenous variables.\n');
        end
        if ~isreal(res) || ~isreal(Y)
            fprintf('\nWARNING: Imaginary parts detected in the residuals or endogenous variables.\n');
        end
        skipline()
    end

    if verbose
        str = sprintf('Iter: %s,\t err. = %s, \t time = %s',num2str(iter),num2str(err), num2str(etime(clock,h2)));
        disp(str);
    end


    if err < options_.dynatol.f
        stop = 1 ;
        break
    end

    iA = iA(1:m,:);
    A = sparse(iA(:,1),iA(:,2),iA(:,3),periods*ny,periods*ny);

    if endogenous_terminal_period && iter>1
        dy = ZERO;
        dy(1:i_rows(end)) = -lin_solve( A(1:i_rows(end),1:i_rows(end)), res(1:i_rows(end)) );
    else
        dy = -lin_solve( A, res );
    end

    Y(i_upd) =   Y(i_upd) + dy;

end

if endogenous_terminal_period
    err = evaluate_max_dynamic_residual(model_dynamic, Y, oo_.exo_simul, params, steady_state, o_periods, ny, max_lag, lead_lag_incidence);
    periods = o_periods;
end


if stop
    if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y)) || ~isreal(res) || ~isreal(Y)
        oo_.deterministic_simulation.status = false;% NaN or Inf occurred
        oo_.deterministic_simulation.error = err;
        oo_.deterministic_simulation.iterations = iter;
        oo_.deterministic_simulation.periods = vperiods(1:iter);
        oo_.endo_simul = reshape(Y,ny,periods+maximum_lag+M_.maximum_lead);
        if verbose
            skipline()
            disp(sprintf('Total time of simulation: %s.', num2str(etime(clock,h1))))
            if ~isreal(res) || ~isreal(Y)
                disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
            else
                disp('Simulation terminated with NaN or Inf in the residuals or endogenous variables.')
            end
            disp('There is most likely something wrong with your model. Try model_diagnostics or another simulation method.')
            printline(105)
        end
    else
        if verbose
            skipline();
            disp(sprintf('Total time of simulation: %s', num2str(etime(clock,h1))))
            printline(56)
        end
        oo_.deterministic_simulation.status = true;% Convergency obtained.
        oo_.deterministic_simulation.error = err;
        oo_.deterministic_simulation.iterations = iter;
        oo_.deterministic_simulation.periods = vperiods(1:iter);
        oo_.endo_simul = reshape(Y,ny,periods+maximum_lag+M_.maximum_lead);
    end
elseif ~stop
    if verbose
        skipline();
        disp(sprintf('Total time of simulation: %s.', num2str(etime(clock,h1))))
        disp('Maximum number of iterations is reached (modify option maxit).')
        printline(62)
    end
    oo_.deterministic_simulation.status = false;% more iterations are needed.
    oo_.deterministic_simulation.error = err;
    oo_.deterministic_simulation.periods = vperiods(1:iter);
    oo_.deterministic_simulation.iterations = options_.simul.maxit;
end

if verbose
    skipline();
end
