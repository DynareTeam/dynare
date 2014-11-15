function [x,fval,exitflag] = simplex_optimization_routine(objective_function,x,options,var_names,varargin)

% Nelder-Mead like optimization routine (see http://en.wikipedia.org/wiki/Nelder-Mead_method)
%
% By default the standard values for the reflection, the expansion, the contraction
% and the shrink coefficients are used (alpha = 1, chi = 2, psi = 1 / 2 and σ = 1 / 2).
%
% The routine automatically restarts from the current solution while amelioration is possible.
%
% INPUTS 
%  o objective_function     [string]                  Name of the objective function to be minimized.
%  o x                      [double]                  n*1 vector, starting guess of the optimization routine.
%  o options                [structure]               Options of this implementation of the simplex algorithm.
%  o var_names              [cell]                    Names of parameters
%                                                       for verbose output
%  o varargin               [cell of structures]      Structures to be passed to the objective function.
%
%     varargin{1} --> DynareDataset
%     varargin{2} --> DatasetInfo
%     varargin{3} --> DynareOptions
%     varargin{4} --> Model
%     varargin{5} --> EstimatedParameters
%     varargin{6} --> BayesInfo
%     varargin{1} --> DynareResults
%
% OUTPUTS 
%  o x                      [double]                  n*1 vector, estimate of the optimal inputs.
%  o fval                   [double]                  scalar, value of the objective at the optimum.
%  o exitflag               [integer]                 scalar equal to 0 or 1 (0 if the algorithm did not converge to
%                                                     a minimum).

% Copyright (C) 2010-2013 Dynare Team
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

% Set verbose mode
verbose = options.verbosity;

% Set number of control variables.
number_of_variables = length(x);

% get options.
if isempty(options.maxfcall)
    max_func_calls = options.maxfcallfactor*number_of_variables
end

% Set tolerance parameter.
if isfield(options,'tolerance') && isfield(options.tolerance,'x')
    x_tolerance = options.tolerance.x;
else
    x_tolerance = 1e-4;
end

% Set tolerance parameter.
if isfield(options,'tolerance') && isfield(options.tolerance,'f')
    f_tolerance = options.tolerance.f;
else
    f_tolerance = 1e-4;
end

% Set maximum number of iterations.
if isfield(options,'maxiter')
    max_iterations = options.maxiter;
else
    max_iterations = 5000;
end

% Set reflection parameter.
if isfield(options,'reflection_parameter')
    if isfield(options.reflection_parameter,'value')
        rho = options.reflection_parameter.value;
    else
        rho = 1.0;
    end
    if isfield(options.reflection_parameter,'random')
        randomize_rho = options.reflection_parameter.random;
        lambda_rho = 1/rho;
    else
        randomize_rho = 0;
    end
else
    rho = 1.0;
    randomize_rho = 0;
end

% Set expansion parameter.
if isfield(options,'expansion_parameter')
    if isfield(options.expansion_parameter,'value')
        chi = options.expansion_parameter.value;
    else
        chi = 2.0;
    end
    if isfield(options.expansion_parameter,'random')
        randomize_chi = options.expansion_parameter.random;
        lambda_chi = 1/chi;
    else
        randomize_chi = 0;
    end
    if isfield(options.expansion_parameter,'optim')
        optimize_expansion_parameter = options.expansion_parameter.optim;
    else
        optimize_expansion_parameter = 0;
    end
else
    chi = 2.0;
    randomize_chi = 0;
    optimize_expansion_parameter = 1;
end

% Set contraction parameter.
if isfield(options,'contraction_parameter')
    if isfield(options.contraction_parameter,'value')
        psi = options.contraction_parameter.value;
    else
        psi = 0.5;
    end
    if isfield(options.contraction_parameter,'random')
        randomize_psi = options.expansion_parameter.random;
    else
        randomize_psi = 0;
    end
else
    psi = 0.5;
    randomize_psi = 0;
end

% Set shrink parameter.
if isfield(options,'shrink_parameter')
    if isfield(options.shrink_parameter,'value')
        sigma = options.shrink_parameter.value;
    else
        sigma = 0.5;
    end
    if isfield(options.shrink_parameter,'random')
        randomize_sigma = options.shrink_parameter.random;
    else
        randomize_sigma = 0;
    end
else
    sigma = 0.5;
    randomize_sigma = 0;
end

% Set delta parameter.
if isfield(options,'delta_parameter')% Size of the simplex
    delta = options.delta_parameter;
else
    delta = 0.05;
end
DELTA = delta;
zero_delta = delta/200;% To be used instead of delta if x(i) is zero.

% Set max_no_improvements.
if isfield(options,'max_no_improvements')
    max_no_improvements = options.max_no_improvements;
else
    max_no_improvements = number_of_variables*10;
end

% Set vector of indices.
unit_vector = ones(1,number_of_variables);
trend_vector_1 = 1:number_of_variables;
trend_vector_2 = 2:(number_of_variables+1);

% Set initial simplex around the initial guess (x).
if verbose
    skipline(3)
    disp('+----------------------+')
    disp(' SIMPLEX INITIALIZATION ')
    disp('+----------------------+')
    skipline()
end
initial_point = x;
[initial_score,junk1,junk2,nopenalty] = feval(objective_function,x,varargin{:});
if ~nopenalty
    error('simplex_optimization_routine:: Initial condition is wrong!')
else
    [v,fv,delta] = simplex_initialization(objective_function,initial_point,initial_score,delta,zero_delta,1,varargin{:});
    func_count = number_of_variables + 1;
    iter_count = 1;
    if verbose
        disp(['Objective function value: ' num2str(fv(1))])
        disp(['Current parameter values: '])
        fprintf(1,'%s: \t\t\t %s \t\t\t %s \t\t\t %s \t\t\t %s \t\t\t %s \n','Names','Best point', 'Worst point', 'Mean values', 'Min values', 'Max values');
        for i=1:number_of_variables
            fprintf(1,'%s: \t\t\t %+8.6f \t\t\t %+8.6f \t\t\t %+8.6f \t\t\t %+8.6f \t\t\t %+8.6f \n',var_names{i},v(i,1), v(i,end), mean(v(i,:),2), min(v(i,:),[],2), max(v(i,:),[],2));
        end
        skipline()
    end
end

vold = v;

no_improvements = 0;
simplex_init = 1;
simplex_iterations = 1;
max_simplex_algo_iterations = 3;
simplex_algo_iterations = 1;
best_point = v(:,1);
best_point_score = fv(1);

convergence = 0;
tooslow = 0;

iter_no_improvement_break = 0;
max_no_improvement_break = 1;

while (func_count < max_func_calls) && (iter_count < max_iterations) && (simplex_algo_iterations<=max_simplex_algo_iterations)
    % Do we really need to continue?
    critF = max(abs(fv(1)-fv(trend_vector_2)));
    critX = max(max(abs(v(:,trend_vector_2)-v(:,unit_vector))));
    if critF <= max(f_tolerance,10*eps(fv(1))) && critX <= max(x_tolerance,10*eps(max(v(:,1))))
        convergence = 1;
    end
    if critX <= x_tolerance^2 && critF>1
        tooslow = 1;
    end
    % Set random reflection and expansion parameters if needed.
    if randomize_rho
        rho = -log(rand)/lambda_rho;
    end
    if randomize_chi
        chi = -log(rand)/lambda_chi;
    end
    % Set random contraction and shrink parameters if needed.
    if randomize_psi
        psi = rand;
    end
    if randomize_sigma
        sigma = rand;
    end
    % Compute the reflection point
    xbar = mean(v(:,trend_vector_1),2); % Average of the n best points.
    xr = xbar + rho*(xbar-v(:,end));
    x  = xr;
    fxr = feval(objective_function,x,varargin{:});
    func_count = func_count+1;
    if fxr < fv(1)% xr is better than previous best point v(:,1).
        % Calculate the expansion point
        xe = xbar + rho*chi*(xbar-v(:,end));
        x  = xe;
        fxe = feval(objective_function,x,varargin{:});
        func_count = func_count+1;
        if fxe < fxr% xe is even better than xr.
            if optimize_expansion_parameter
                if verbose>1
                    skipline(2)
                    disp('Compute optimal expansion...')
                end
                xee  = xbar + rho*chi*1.01*(xbar-v(:,end));
                x    = xee;
                fxee = feval(objective_function,x,varargin{:});
                func_count = func_count+1;
                if fxee<fxe
                    decrease = 1;
                    weight = rho*chi*1.02;
                    fxeee_old = fxee;
                    xeee_old  = xee;
                    if verbose>1
                        fprintf(1,'Weight =        ');
                    end
                    while decrease
                        weight = 1.02*weight;
                        if verbose>1
                            fprintf(1,'\b\b\b\b\b\b\b %6.4f',weight);
                        end
                        xeee   = xbar + weight*(xbar-v(:,end));
                        x      = xeee;
                        fxeee  = feval(objective_function,x,varargin{:});
                        func_count = func_count+1;
                        if (fxeee<fxeee_old) && -(fxeee-fxeee_old)>f_tolerance*10*fxeee_old
                            fxeee_old = fxeee;
                            xeee_old  = xeee;
                        else
                            decrease = 0;
                        end
                    end
                    if verbose>1
                        fprintf(1,'\n');
                    end
                    xe  = xeee_old;
                    fxe = fxeee_old;
                else
                    decrease = 1;
                    weight = rho*chi;
                    fxeee_old = fxee;
                    xeee_old  = xee;
                    if verbose>1
                        fprintf(1,'Weight =        ');
                    end
                    while decrease
                        weight = weight/1.02;
                        if verbose>1
                            fprintf(1,'\b\b\b\b\b\b\b %6.4f',weight);
                        end
                        xeee   = xbar + weight*(xbar-v(:,end));
                        x      = xeee;
                        fxeee  = feval(objective_function,x,varargin{:});
                        func_count = func_count+1;
                        if (fxeee<fxeee_old) && -(fxeee-fxeee_old)>f_tolerance*10*fxeee_old
                            fxeee_old = fxeee;
                            xeee_old  = xeee;
                        else
                            decrease = 0;
                        end
                    end
                    if verbose>1
                        fprintf(1,'\n');
                    end
                    xe  = xeee_old;
                    fxe = fxeee_old;
                end
                if verbose>1
                    disp('Done!')
                    skipline(2)
                end
            end
            v(:,end) = xe;
            fv(end)  = fxe;
            move = 'expand';
        else% if xe is not better than xr.
            v(:,end) = xr;
            fv(end)  = fxr;
            move = 'reflect-1';
        end
    else% xr is not better than previous best point v(:,1).
        if fxr < fv(number_of_variables)% xr is better than previous point v(:,n).
            v(:,end) = xr;
            fv(end) = fxr;
            move = 'reflect-0';
        else% xr is not better than previous point v(:,n).
            if fxr < fv(end)% xr is better than previous worst point [=> outside contraction].
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                x  = xc;
                fxc = feval(objective_function,x,varargin{:});
                func_count = func_count+1;
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(end) = fxc;
                    move = 'contract outside';
                else
                    move = 'shrink';
                end
            else% xr is the worst point [=> inside contraction].
                xcc = (1-psi)*xbar + psi*v(:,end);
                x   = xcc;
                fxcc = feval(objective_function,x,varargin{:});
                func_count = func_count+1;
                if fxcc < fv(end)
                    v(:,end) = xcc;
                    fv(end) = fxcc;
                    move = 'contract inside';
                else
                    % perform a shrink
                    move = 'shrink';
                end
            end
            if strcmp(move,'shrink')
                for j=trend_vector_2
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    x = v(:,j);
                    fv(j) = feval(objective_function,x,varargin{:});
                end
                func_count = func_count + number_of_variables;
            end
        end
    end
    % Sort n+1 points by incresing order of the objective function values.
    [fv,sort_idx] = sort(fv);
    v = v(:,sort_idx);
    iter_count = iter_count + 1;
    simplex_iterations = simplex_iterations+1;
    if verbose>1
        disp(['Simplex iteration number: ' int2str(simplex_iterations) '-' int2str(simplex_init) '-' int2str(simplex_algo_iterations)])
        disp(['Simplex move:             ' move])
        disp(['Objective function value: ' num2str(fv(1))])
        disp(['Mode improvement:         ' num2str(best_point_score-fv(1))])
        disp(['Norm of dx:               ' num2str(norm(best_point-v(:,1)))])
        disp(['Norm of dSimplex:         ' num2str(norm(v-vold,'fro'))])
        disp(['Crit. f:                  ' num2str(critF)])
        disp(['Crit. x:                  ' num2str(critX)])
        skipline()
    end
    if verbose && max(abs(best_point-v(:,1)))>x_tolerance;
        if verbose<2
            disp(['Simplex iteration number: ' int2str(simplex_iterations) '-' int2str(simplex_init) '-' int2str(simplex_algo_iterations)])
            disp(['Objective function value: ' num2str(fv(1))])
            disp(['Mode improvement:         ' num2str(best_point_score-fv(1))])
            disp(['Norm of dx:               ' num2str(norm(best_point-v(:,1)))])
            disp(['Norm of dSimplex:         ' num2str(norm(v-vold,'fro'))])
            disp(['Crit. f:                  ' num2str(critF)])
            disp(['Crit. x:                  ' num2str(critX)])
            skipline()
        end
        disp(['Current parameter values: '])
        fprintf(1,'%s: \t\t\t %s \t\t\t %s \t\t\t %s \t\t\t %s \t\t\t %s \n','Names','Best point', 'Worst point', 'Mean values', 'Min values', 'Max values');
        for i=1:number_of_variables
            fprintf(1,'%s: \t\t\t %+8.6f \t\t\t %+8.6f \t\t\t %+8.6f \t\t\t %+8.6f \t\t\t %+8.6f \n',var_names{i}, v(i,1), v(i,end), mean(v(i,:),2), min(v(i,:),[],2), max(v(i,:),[],2));
        end
        skipline()
    end
    if abs(best_point_score-fv(1))<f_tolerance
        no_improvements = no_improvements+1;
    else
        no_improvements = 0;
    end
    best_point = v(:,1);
    best_point_score = fv(1);
    vold = v;
    if no_improvements>max_no_improvements
        if verbose
            disp(['NO SIGNIFICANT IMPROVEMENT AFTER ' int2str(no_improvements) ' ITERATIONS!'])
        end
        if simplex_algo_iterations<=max_simplex_algo_iterations
            % Compute the size of the simplex
            delta = delta*1.05;
            % Compute the new initial simplex.
            [v,fv,delta] = simplex_initialization(objective_function,best_point,best_point_score,delta,zero_delta,1,varargin{:});
            if verbose
                disp(['(Re)Start with a lager simplex around the based on the best current '])
                disp(['values for the control variables. '])
                disp(['New value of delta (size of the new simplex) is: '])
                for i=1:number_of_variables
                    fprintf(1,'%s: \t\t\t %+8.6f \n',var_names{i}, delta(i));
                end
            end
            % Reset counters
            no_improvements = 0;
            func_count = func_count + number_of_variables;
            iter_count = iter_count+1;
            iter_no_improvement_break = iter_no_improvement_break + 1;
            simplex_init = simplex_init+1;
            simplex_iterations = simplex_iterations+1;
            skipline(2)
        end
    end
    if ((func_count==max_func_calls) || (iter_count==max_iterations) || (iter_no_improvement_break==max_no_improvement_break) || convergence || tooslow)
        [v,fv,delta] = simplex_initialization(objective_function,best_point,best_point_score,DELTA,zero_delta,1,varargin{:});
        if func_count==max_func_calls
            if verbose
                disp(['MAXIMUM NUMBER OF OBJECTIVE FUNCTION CALLS EXCEEDED (' int2str(max_func_calls) ')!'])
            end
        elseif iter_count== max_iterations
            if verbose
                disp(['MAXIMUM NUMBER OF ITERATIONS EXCEEDED (' int2str(max_iterations) ')!'])
            end
        elseif iter_no_improvement_break==max_no_improvement_break
            if verbose
                disp(['MAXIMUM NUMBER OF SIMPLEX REINITIALIZATION EXCEEDED (' int2str(max_no_improvement_break) ')!'])
            end
            iter_no_improvement_break = 0;
            if simplex_algo_iterations==max_simplex_algo_iterations
                max_no_improvements = Inf;% Do not stop until convergence is reached!
                continue
            end
        elseif tooslow
            disp(['CONVERGENCE NOT ACHIEVED AFTER ' int2str(simplex_iterations) ' ITERATIONS! IMPROVING TOO SLOWLY!'])
        else
            disp(['CONVERGENCE ACHIEVED AFTER ' int2str(simplex_iterations) ' ITERATIONS!'])
        end
        if simplex_algo_iterations<max_simplex_algo_iterations
            % Compute the size of the simplex
            delta = delta*1.05;
            % Compute the new initial simplex.
            [v,fv,delta] = simplex_initialization(objective_function,best_point,best_point_score,delta,zero_delta,1,varargin{:});
            if verbose
                disp(['(Re)Start with a lager simplex around the based on the best current '])
                disp(['values for the control variables. '])
                disp(['New value of delta (size of the new simplex) is: '])
                for i=1:number_of_variables
                    fprintf(1,'%s: \t\t\t %+8.6f \n',var_names{i}, delta(i));
                end
            end
            % Reset counters
            func_count=0;
            iter_count=0;
            convergence = 0;
            no_improvements = 0;
            func_count = func_count + number_of_variables;
            iter_count = iter_count+1;
            simplex_iterations = simplex_iterations+1;
            simplex_algo_iterations = simplex_algo_iterations+1;
            skipline(2)
        else
            break
        end
   end
end% while loop.

x(:) = v(:,1);
fval = fv(1);
exitflag = 1;

if func_count>= max_func_calls
    disp('Maximum number of objective function calls has been exceeded!')
    exitflag = 0;
end

if iter_count>= max_iterations
    disp('Maximum number of iterations has been exceeded!')
    exitflag = 0;
end




function [v,fv,delta] = simplex_initialization(objective_function,point,point_score,delta,zero_delta,check_delta,varargin)
    n = length(point);
    v  = zeros(n,n+1);
    v(:,1) = point;
    fv = zeros(n+1,1);
    fv(1) = point_score;
    if length(delta)==1
        delta = repmat(delta,n,1);
    end
    for j = 1:n
        y = point;
        if y(j) ~= 0
            y(j) = (1 + delta(j))*y(j);
        else
            y(j) = zero_delta;
        end
        v(:,j+1) = y;
        x = y;
        [fv(j+1),junk1,junk2,nopenalty_flag] = feval(objective_function,x,varargin{:});
        if check_delta
            while ~nopenalty_flag
                if y(j)~=0
                    delta(j) = delta(j)/1.1;
                else
                    zero_delta = zero_delta/1.1;
                end
                y = point;
                if y(j) ~= 0
                    y(j) = (1 + delta(j))*y(j);
                else
                    y(j) = zero_delta;
                end
                v(:,j+1) = y;
                x = y;
                [fv(j+1),junk1,junk2,nopenalty_flag] = feval(objective_function,x,varargin{:});
            end
        end
    end
    % Sort by increasing order of the objective function values.
    [fv,sort_idx] = sort(fv);
    v = v(:,sort_idx);