function [endogenousvariables, info] = sim1(endogenousvariables, exogenousvariables, steadystate, M, options)

% Performs deterministic simulations with lead or lag on one period. Uses sparse matrices.
%
% INPUTS
%   - endogenousvariables [double] N*T array, paths for the endogenous variables (initial guess).
%   - exogenousvariables  [double] T*M array, paths for the exogenous variables.
%   - steadystate       [double] N*1 array, steady state for the endogenous variables.
%   - M                   [struct] contains a description of the model.
%   - options             [struct] contains various options.
% OUTPUTS
%   - endogenousvariables [double] N*T array, paths for the endogenous variables (solution of the perfect foresight model).
%   - info                [struct] contains informations about the results.
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 1996-2016 Dynare Team
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

verbose = options.verbosity;

endogenous_terminal_period = options.endogenous_terminal_period;
vperiods = options.periods*ones(1,options.simul.maxit);
azero = options.dynatol.f/1e7;

lead_lag_incidence = M.lead_lag_incidence;

ny = M.endo_nbr;

maximum_lag = M.maximum_lag;
max_lag = M.maximum_endo_lag;

nyp = nnz(lead_lag_incidence(1,:)) ;
ny0 = nnz(lead_lag_incidence(2,:)) ;
nyf = nnz(lead_lag_incidence(3,:)) ;

nd = nyp+ny0+nyf;
stop = 0 ;

periods = options.periods;
params = M.params;

i_cols_1 = nonzeros(lead_lag_incidence(2:3,:)');
i_cols_A1 = find(lead_lag_incidence(2:3,:)');
i_cols_A1 = i_cols_A1(:);
i_cols_T = nonzeros(lead_lag_incidence(1:2,:)');
i_cols_0 = nonzeros(lead_lag_incidence(2,:)');
i_cols_A0 = find(lead_lag_incidence(2,:)');
i_cols_A0 = i_cols_A0(:);
i_cols_j = (1:nd)';
i_upd = maximum_lag*ny+(1:periods*ny);

Y = endogenousvariables(:);

if verbose
    skipline()
    printline(56)
    disp('MODEL SIMULATION:')
    skipline()
end

model_dynamic = str2func([M.fname,'_dynamic']);
z = Y(find(lead_lag_incidence'));

[d1,jacobian] = model_dynamic(z, exogenousvariables, params, steadystate,maximum_lag+1);

res = zeros(periods*ny,1);

o_periods = periods;

if endogenous_terminal_period
    ZERO = zeros(length(i_upd),1);
end

h1 = clock ;
iA = zeros(periods*M.NNZDerivatives(1),3);

for iter = 1:options.simul.maxit
    h2 = clock ;
    i_rows = (1:ny)';
    i_cols_A = find(lead_lag_incidence');
    i_cols_A = i_cols_A(:);
    i_cols = i_cols_A+(maximum_lag-1)*ny;
    m = 0;
    for it = (maximum_lag+1):(maximum_lag+periods)
        [d1,jacobian] = model_dynamic(Y(i_cols), exogenousvariables, params, steadystate,it);
        if it == maximum_lag+periods && it == maximum_lag+1
            [r,c,v] = find(jacobian(:,i_cols_0));
            iA((1:length(v))+m,:) = [i_rows(r(:)),i_cols_A0(c(:)),v(:)];
        elseif it == maximum_lag+periods
            [r,c,v] = find(jacobian(:,i_cols_T));
            iA((1:length(v))+m,:) = [i_rows(r(:)),i_cols_A(i_cols_T(c(:))),v(:)];
        elseif it == maximum_lag+1
            [r,c,v] = find(jacobian(:,i_cols_1));
            iA((1:length(v))+m,:) = [i_rows(r(:)),i_cols_A1(c(:)),v(:)];
        else
            [r,c,v] = find(jacobian(:,i_cols_j));
            iA((1:length(v))+m,:) = [i_rows(r(:)),i_cols_A(c(:)),v(:)];
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
    if options.debug
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
    if err < options.dynatol.f
        stop = 1 ;
        break
    end
    iA = iA(1:m,:);
    A = sparse(iA(:,1),iA(:,2),iA(:,3),periods*ny,periods*ny);

    if iter>1,
        dy0=dy; % store previous update step
    end
    
    if endogenous_terminal_period && iter>1
        dy = ZERO;
        if options.simul.robust_lin_solve
            dy(1:i_rows(end)) = -lin_solve_robust( A(1:i_rows(end),1:i_rows(end)), res(1:i_rows(end)),verbose );            
        else
            dy(1:i_rows(end)) = -lin_solve( A(1:i_rows(end),1:i_rows(end)), res(1:i_rows(end)), verbose );
        end
    else
        if options.simul.robust_lin_solve
            dy = -lin_solve_robust( A, res, verbose );            
        else
            dy = -lin_solve( A, res, verbose );
    end
    if iter>1 && (any(~isreal(dy)) || any(isnan(dy)) || any(isinf(dy)) || (err/err0)>10),
        % relaxation step, reducing length of update of previous
        % iteration
        if verbose
            display_critical_variables(reshape(dy,[ny periods])', M);
        end
        disp('Try relaxation step.')
        Y=Y0;
        dy=dy0/10; %max(10,(err/err0));
        err=err0;
    end
    Y0=Y; % store previous endogenous variables path
    err0=err;
    
    Y(i_upd) = Y(i_upd) + dy;
end

if endogenous_terminal_period
    err = evaluate_max_dynamic_residual(model_dynamic, Y, exogenousvariables, params, steadystate, o_periods, ny, max_lag, lead_lag_incidence);
    periods = o_periods;
end


if stop
    if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y)) || ~isreal(res) || ~isreal(Y)
        info.status = false;% NaN or Inf occurred
        info.error = err;
        info.iterations = iter;
        info.periods = vperiods(1:iter);
        endogenousvariables = reshape(Y,ny,periods+maximum_lag+M.maximum_lead);
        if verbose
            skipline()
            disp(sprintf('Total time of simulation: %s.', num2str(etime(clock,h1))))
            if ~isreal(res) || ~isreal(Y)
                disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
            else
                disp('Simulation terminated with NaN or Inf in the residuals or endogenous variables.')
            end
            display_critical_variables(reshape(dy,[ny periods])', M);
            disp('There is most likely something wrong with your model. Try model_diagnostics or another simulation method.')
            printline(105)
        end
    else
        if verbose
            skipline();
            disp(sprintf('Total time of simulation: %s', num2str(etime(clock,h1))))
            printline(56)
        end
        info.status = true;% Convergency obtained.
        info.error = err;
        info.iterations = iter;
        info.periods = vperiods(1:iter);
        endogenousvariables = reshape(Y,ny,periods+maximum_lag+M.maximum_lead);
    end
elseif ~stop
    if verbose
        skipline();
        disp(sprintf('Total time of simulation: %s.', num2str(etime(clock,h1))))
        disp('Maximum number of iterations is reached (modify option maxit).')
        printline(62)
    end
    info.status = false;% more iterations are needed.
    info.error = err;
    info.periods = vperiods(1:iter);
    info.iterations = options.simul.maxit;
end

if verbose
    skipline();
end

function x = lin_solve( A, b,verbose)
    if norm( b ) < sqrt( eps ) % then x = 0 is a solution
        x = 0;
        return
    end
    
    x = A\b;
    x( ~isfinite( x ) ) = 0;
    relres = norm( b - A * x ) / norm( b );
    if relres > 1e-6 && verbose
        fprintf( 'WARNING : Failed to find a solution to the linear system.\n' );
    end
    
function [ x, flag, relres ] = lin_solve_robust( A, b , verbose)
    if norm( b ) < sqrt( eps ) % then x = 0 is a solution
        x = 0;
        flag = 0;
        relres = 0;
        return
    end
    
    x = A\b;
    x( ~isfinite( x ) ) = 0;
    [ x, flag, relres ] = bicgstab( A, b, [], [], [], [], x ); % returns immediately if x is a solution
    if flag == 0
        return
    end

    disp( relres );

    if verbose
        fprintf( 'Initial bicgstab failed, trying alternative start point.\n' );
    end
    old_x = x;
    old_relres = relres;
    [ x, flag, relres ] = bicgstab( A, b );
    if flag == 0
        return
    end

    if verbose
        fprintf( 'Alternative start point also failed with bicgstab, trying gmres.\n' );
    end
    if old_relres < relres
        x = old_x;
    end
    [ x, flag, relres ] = gmres( A, b, [], [], [], [], [], x );
    if flag == 0
        return
    end

    if verbose
        fprintf( 'Initial gmres failed, trying alternative start point.\n' );
    end
    old_x = x;
    old_relres = relres;
    [ x, flag, relres ] = gmres( A, b );
    if flag == 0
        return
    end

    if verbose
        fprintf( 'Alternative start point also failed with gmres, using the (SLOW) Moore-Penrose Pseudo-Inverse.\n' );
    end
    if old_relres < relres
        x = old_x;
        relres = old_relres;
    end
    old_x = x;
    old_relres = relres;
    x = pinv( full( A ) ) * b;
    relres = norm( b - A * x ) / norm( b );
    if old_relres < relres
        x = old_x;
        relres = old_relres;
    end
    flag = relres > 1e-6;
    if flag ~= 0 && verbose
        fprintf( 'WARNING : Failed to find a solution to the linear system\n' );
    end
    
function display_critical_variables(dyy, M)

            if any(isnan(dyy))
                indx = find(any(isnan(dyy)));
                endo_names=cellstr(M.endo_names(indx,:));
                disp('Last iteration provided NaN for the following variables:')
                fprintf('%s, ',endo_names{:}),
                fprintf('\n'),
            end
            if any(isinf(dyy))
                indx = find(any(isinf(dyy)));
                endo_names=cellstr(M.endo_names(indx,:));
                disp('Last iteration diverged (Inf) for the following variables:')
                fprintf('%s, ',endo_names{:}),
                fprintf('\n'),
            end
            if any(~isreal(dyy))
                indx = find(any(~isreal(dyy)));
                endo_names=cellstr(M.endo_names(indx,:));
                disp('Last iteration provided complex number for the following variables:')
                fprintf('%s, ',endo_names{:}),
                fprintf('\n'),                
            end

        