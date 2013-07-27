function simul()
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

test_for_deep_parameters_calibration(M_);

if options_.stack_solve_algo < 0 || options_.stack_solve_algo > 6
    error('SIMUL: stack_solve_algo must be between 0 and 6')
end

if ~options_.block && ~options_.bytecode && options_.stack_solve_algo ~= 0 ...
        && options_.stack_solve_algo ~= 6
    error('SIMUL: you must use stack_solve_algo=0 or stack_solve_algo=6 when not using block nor bytecode option')
end

if options_.block && ~options_.bytecode && options_.stack_solve_algo == 5
    error('SIMUL: you can''t use stack_solve_algo = 5 without bytecode option')
end

if (options_.block || options_.bytecode) && options_.stack_solve_algo == 6
    error('SIMUL: you can''t use stack_solve_algo = 6 with block or bytecode option')
end

if exist('OCTAVE_VERSION') && options_.stack_solve_algo == 2
    error('SIMUL: you can''t use stack_solve_algo = 2 under Octave')
end

if size(M_.lead_lag_incidence,2)-nnz(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)) > 0
    mess = ['SIMUL: error in model specification : variable ' M_.endo_names(find(M_.lead_lag_incidence(M_.maximum_lag+1,:)==0),:)] ;
    mess = [mess ' doesn''t appear as current variable.'] ; 
    error (mess) ;
end

if options_.periods == 0
    error('SIMUL: number of periods for the simulation isn''t specified')
end

if ~ options_.initval_file
    if isempty(options_.datafile)
        make_ex_;
        make_y_;
    else
        read_data_;
    end
end

if isempty(options_.scalv) || options_.scalv == 0
    options_.scalv = oo_.steady_state ;
end

options_.scalv= 1 ;

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

if(options_.block)
    if(options_.bytecode)
        [info, oo_.endo_simul] = bytecode('dynamic');
        if info == 1
            oo_.deterministic_simulation.status = 0;
        else
            oo_.deterministic_simulation.status = 1;
        end;
        mexErrCheck('bytecode', info);
    else
        eval([M_.fname '_dynamic']);
    end;
else
    if(options_.bytecode)
        [info, oo_.endo_simul]=bytecode('dynamic');
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
    end;
end;

dyn2vec;