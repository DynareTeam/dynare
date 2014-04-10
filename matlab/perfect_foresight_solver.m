function perfect_foresight_solver()
% Computes deterministic simulations
%  
% INPUTS
%   None
%  
% OUTPUTS
%   none
%    
% ALGORITHM
%   
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 1996-2014 Dynare Team
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

if options_.stack_solve_algo < 0 || options_.stack_solve_algo > 6
    error('PERFECT_FORESIGHT_SOLVER: stack_solve_algo must be between 0 and 6')
end

if ~options_.block && ~options_.bytecode && options_.stack_solve_algo ~= 0 ...
        && options_.stack_solve_algo ~= 6
    error('PERFECT_FORESIGHT_SOLVER: you must use stack_solve_algo=0 or stack_solve_algo=6 when not using block nor bytecode option')
end

if options_.block && ~options_.bytecode && options_.stack_solve_algo == 5
    error('PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 5 without bytecode option')
end

if (options_.block || options_.bytecode) && options_.stack_solve_algo == 6
    error('PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 6 with block or bytecode option')
end

if isoctave && options_.stack_solve_algo == 2
    error('PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 2 under Octave')
end


if isempty(oo_.endo_simul) || any(size(oo_.endo_simul) ~= [ M_.endo_nbr, M_.maximum_endo_lag+options_.periods+M_.maximum_endo_lead ])
    error('PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
end

if isempty(oo_.exo_simul) || any(size(oo_.exo_simul) ~= [ M_.maximum_lag+options_.periods+M_.maximum_lead, M_.exo_nbr ])
    error('PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
end


if isempty(options_.scalv) || options_.scalv == 0
    options_.scalv = oo_.steady_state;
end

options_.scalv= 1;

if options_.debug
    model_static = str2func([M_.fname,'_static']);
    for ii=1:size(oo_.exo_simul,1)
        [residual(:,ii)] = model_static(oo_.steady_state, oo_.exo_simul(ii,:),M_.params);
    end
    problematic_periods=find(any(isinf(residual)) | any(isnan(residual)))-M_.maximum_endo_lag;
    if ~isempty(problematic_periods) 
        period_string=num2str(problematic_periods(1));
        for ii=2:length(problematic_periods)
            period_string=[period_string, ', ', num2str(problematic_periods(ii))];
        end
        fprintf('\n\nWARNING: Value for the exogenous variable(s) in period(s) %s inconsistent with the static model.\n',period_string);   
        fprintf('WARNING: Check for division by 0.\n')
    end
end

% Effectively compute simulation, possibly with homotopy
if options_.no_homotopy
    simulation_core;
else
    exosim = oo_.exo_simul;
    exoinit = repmat(oo_.exo_steady_state',M_.maximum_lag+options_.periods+M_.maximum_lead,1);
    endosim = oo_.endo_simul;
    endoinit = repmat(oo_.steady_state, 1,M_.maximum_endo_lag+options_.periods+M_.maximum_endo_lead);

    current_weight = 0; % Current weight of target point in convex combination
    step = 1;
    success_counter = 0;

    while (step > options_.dynatol.x)

        new_weight = current_weight + step; % Try this weight, and see if it succeeds
        if new_weight >= 1
            new_weight = 1; % Don't go beyond target point
            step = new_weight - current_weight;
        end

        % Compute convex combination for exo path and initial/terminal endo conditions
        % But take care of not overwriting the computed part of oo_.endo_simul
        oo_.exo_simul = exosim*new_weight + exoinit*(1-new_weight);
        endocombi = endosim*new_weight + endoinit*(1-new_weight);
        oo_.endo_simul(:,1:M_.maximum_endo_lag) = endocombi(:,1:M_.maximum_endo_lag);
        oo_.endo_simul(:,(end-M_.maximum_endo_lead):end) = endocombi(:,(end-M_.maximum_endo_lead):end);

        saved_endo_simul = oo_.endo_simul;

        simulation_core;

        if oo_.deterministic_simulation.status == 1
            current_weight = new_weight;
            if current_weight >= 1
                break
            end
            success_counter = success_counter + 1;
            if success_counter >= 3
                success_counter = 0;
                step = step * 2;
                disp([ 'Homotopy step succeeded, doubling step size (completed ' sprintf('%.1f', current_weight*100) '%, step size ' sprintf('%.3g', step) ')' ])
            else
                disp([ 'Homotopy step succeeded (completed ' sprintf('%.1f', current_weight*100) '%, step size ' sprintf('%.3g', step) ')' ])
            end
        else
            oo_.endo_simul = saved_endo_simul;
            success_counter = 0;
            step = step / 2;
            disp([ 'Homotopy step failed, halving step size (completed ' sprintf('%.1f', current_weight*100) '%, step size ' sprintf('%.3g', step) ')' ])
        end
    end
end

if oo_.deterministic_simulation.status == 1
    disp('Perfect foresight solution found.')
else
    disp('Failed to solve perfect foresight model')
end

dyn2vec;

ts = dseries(transpose(oo_.endo_simul),options_.initial_period,cellstr(M_.endo_names));
assignin('base', 'Simulated_time_series', ts);

end


function simulation_core()

global M_ oo_ options_

if(options_.block)
    if(options_.bytecode)
        [info, oo_.endo_simul] = bytecode('dynamic');
        if info == 1
            oo_.deterministic_simulation.status = 0;
        else
            oo_.deterministic_simulation.status = 1;
        end
        mexErrCheck('bytecode', info);
    else
        eval([M_.fname '_dynamic']);
    end
else
    if(options_.bytecode)
        [info, oo_.endo_simul]=bytecode('dynamic');
        if info == 1
            oo_.deterministic_simulation.status = 0;
        else
            oo_.deterministic_simulation.status = 1;
        end;
        mexErrCheck('bytecode', info);
    else
        if M_.maximum_endo_lead == 0 % Purely backward model
            sim1_purely_backward;
        elseif M_.maximum_endo_lag == 0 % Purely forward model
            sim1_purely_forward;
        else % General case
            if options_.stack_solve_algo == 0
                sim1;
            else % stack_solve_algo = 6
                sim1_lbj;
            end
        end
    end
end

end
