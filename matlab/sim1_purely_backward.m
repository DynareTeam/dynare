function sim1_purely_backward()
% Performs deterministic simulation of a purely backward model

% Copyright (C) 2012-2014 Dynare Team
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
    if size(M_.lead_lag_incidence,1) > 1
        ny0 = nnz(M_.lead_lag_incidence(2,:)); % Number of variables at current period
        nyb = nnz(M_.lead_lag_incidence(1,:)); % Number of variables at previous period
        iyb = find(M_.lead_lag_incidence(1,:)>0); % Indices of variables at previous period
    else
        ny0 = nnz(M_.lead_lag_incidence(1,:)); % Number of variables at current period
        nyb = 0;
        iyb = [];
    end
        

    if ny0 ~= M_.endo_nbr
        error('SIMUL: all endogenous variables must appear at the current period')
    end
    
    model_dynamic = str2func([M_.fname,'_dynamic']);

    oo_.deterministic_simulation.status = 1;

    for it = 2:options_.periods+1
        yb = oo_.endo_simul(:,it-1); % Values at previous period, also used as guess value for current period
        yb1 = yb(iyb);
       
        [tmp, check] = solve1(model_dynamic, [yb1; yb], 1:M_.endo_nbr, nyb+1:nyb+ ...
                              M_.endo_nbr, 1, options_.gstep, ...
                              options_.solve_tolf,options_.solve_tolx, ...
                              options_.simul.maxit,options_.debug,oo_.exo_simul, ...
                              M_.params, oo_.steady_state, it);

        if check
            oo_.deterministic_simulation.status = 0;
        end

        oo_.endo_simul(:,it) = tmp(nyb+1:nyb+M_.endo_nbr);
    end
    