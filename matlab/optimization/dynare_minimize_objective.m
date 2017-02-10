function [opt_par_values,fval,exitflag,hessian_mat,options_,Scale,new_rat_hess_info]=dynare_minimize_objective(objective_function,start_par_value,minimizer_algorithm,options_,bounds,parameter_names,prior_information,Initial_Hessian,varargin)
% function [opt_par_values,fval,exitflag,hessian_mat,options_,Scale,new_rat_hess_info]=dynare_minimize_objective(objective_function,start_par_value,minimizer_algorithm,options_,bounds,parameter_names,prior_information,Initial_Hessian,new_rat_hess_info,varargin)
% Calls a minimizer
%
% INPUTS
%   objective_function  [function handle]                   handle to the objective function
%   start_par_value     [n_params by 1] vector of doubles   starting values for the parameters
%   minimizer_algorithm [scalar double, or string]          code of the optimizer algorithm, or string for the name of a user defined optimization routine (not shipped with dynare).
%   options_            [matlab structure]                  Dynare options structure
%   bounds              [n_params by 2] vector of doubles   2 row vectors containing lower and upper bound for parameters
%   parameter_names     [n_params by 1] cell array          strings containing the parameters names   
%   prior_information   [matlab structure]                  Dynare prior information structure (bayestopt_) provided for algorithm 6
%   Initial_Hessian     [n_params by n_params] matrix       initial hessian matrix provided for algorithm 6
%   new_rat_hess_info   [matlab structure]                  step size info used by algorith 5
%   varargin            [cell array]                        Input arguments for objective function
%    
% OUTPUTS
%   opt_par_values      [n_params by 1] vector of doubles   optimal parameter values minimizing the objective
%   fval                [scalar double]                     value of the objective function at the minimum
%   exitflag            [scalar double]                     return code of the respective optimizer
%   hessian_mat         [n_params by n_params] matrix       hessian matrix at the mode returned by optimizer
%   options_            [matlab structure]                  Dynare options structure (to return options set by algorithms 5)
%   Scale               [scalar double]                     scaling parameter returned by algorith 6
%    
% SPECIAL REQUIREMENTS
%   none.
%  
% 
% Copyright (C) 2014-2016 Dynare Team
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


%% set bounds and parameter names if not already set
n_params=size(start_par_value,1);
if isempty(bounds)
    bounds=[-Inf(n_params,1) Inf(n_params,1)];
end

if isempty(parameter_names)
    parameter_names=[repmat('parameter ',n_params,1),num2str((1:n_params)')];
end

%% initialize function outputs
hessian_mat=[];
Scale=[];
exitflag=1;
fval=NaN;
opt_par_values=NaN(size(start_par_value));
new_rat_hess_info=[];

switch minimizer_algorithm
  case 1
    if isoctave
        error('Optimization algorithm 1 is not available under Octave')
    elseif ~user_has_matlab_license('optimization_toolbox')
        error('Optimization algorithm 1 requires the Optimization Toolbox')
    end
    % Set default optimization options for fmincon.
    optim_options = optimset('display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-6);
    if ~isempty(options_.optim_opt)
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if options_.silent_optimizer
        optim_options = optimset(optim_options,'display','off');
    end
    if options_.analytic_derivation,
        optim_options = optimset(optim_options,'GradObj','on','TolX',1e-7);
    end
    [opt_par_values,fval,exitflag,output,lamdba,grad,hessian_mat] = ...
        fmincon(objective_function,start_par_value,[],[],[],[],bounds(:,1),bounds(:,2),[],optim_options,varargin{:});
  case 2
    %simulating annealing
    sa_options = options_.saopt;
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'neps'
                sa_options.neps = options_list{i,2};
              case 'rt'
                sa_options.rt = options_list{i,2};
              case 'MaxIter'
                sa_options.MaxIter = options_list{i,2};
              case 'TolFun'
                sa_options.TolFun = options_list{i,2};
              case 'verbosity'
                sa_options.verbosity = options_list{i,2};
              case 'initial_temperature'
                sa_options.initial_temperature = options_list{i,2};
              case 'ns'
                sa_options.ns = options_list{i,2};
              case 'nt'
                sa_options.nt = options_list{i,2};
              case 'step_length_c'
                sa_options.step_length_c = options_list{i,2};
              case 'initial_step_length'
                sa_options.initial_step_length = options_list{i,2};
              otherwise
                warning(['solveopt: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    if options_.silent_optimizer
        sa_options.verbosity = 0;
    end
    npar=length(start_par_value);
    [LB, UB]=set_bounds_to_finite_values(bounds, options_.huge_number);
    if sa_options.verbosity
        fprintf('\nNumber of parameters= %d, initial temperatur= %4.3f \n', npar,sa_options.initial_temperature);
        fprintf('rt=  %4.3f; TolFun=  %4.3f; ns=  %4.3f;\n',sa_options.rt,sa_options.TolFun,sa_options.ns);
        fprintf('nt=  %4.3f; neps=  %4.3f; MaxIter=  %d\n',sa_options.nt,sa_options.neps,sa_options.MaxIter);
        fprintf('Initial step length(vm): %4.3f; step_length_c: %4.3f\n', sa_options.initial_step_length,sa_options.step_length_c);
        fprintf('%-20s  %-6s    %-6s    %-6s\n','Name:', 'LB;','Start;','UB;');
        for pariter=1:npar
            fprintf('%-20s  %6.4f;   %6.4f;  %6.4f;\n',parameter_names{pariter}, LB(pariter),start_par_value(pariter),UB(pariter));
        end
    end
    sa_options.initial_step_length= sa_options.initial_step_length*ones(npar,1); %bring step length to correct vector size
    sa_options.step_length_c= sa_options.step_length_c*ones(npar,1); %bring step_length_c to correct vector size
    [opt_par_values, fval,exitflag, n_accepted_draws, n_total_draws, n_out_of_bounds_draws, t, vm] =...
        simulated_annealing(objective_function,start_par_value,sa_options,LB,UB,varargin{:});
  case 3
    if isoctave && ~user_has_octave_forge_package('optim')
        try
            pkg load optim
        catch
            error('Optimization algorithm 3 requires the optim package')
        end
    elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
        error('Optimization algorithm 3 requires the Optimization Toolbox')
    end
    % Set default optimization options for fminunc.
    optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
    if ~isempty(options_.optim_opt)
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if options_.analytic_derivation,
        optim_options = optimset(optim_options,'GradObj','on');
    end
    if options_.silent_optimizer
        optim_options = optimset(optim_options,'display','off');
    end
    if ~isoctave
        [opt_par_values,fval,exitflag] = fminunc(objective_function,start_par_value,optim_options,varargin{:});
    else
        % Under Octave, use a wrapper, since fminunc() does not have a 4th arg
        func = @(x) objective_function(x,varargin{:});
        [opt_par_values,fval,exitflag] = fminunc(func,start_par_value,optim_options);
    end
  case 4
    % Set default options.
    H0 = 1e-4*eye(n_params);
    crit = options_.csminwel.tolerance.f;
    nit = options_.csminwel.maxiter;
    numgrad = options_.gradient_method;
    epsilon = options_.gradient_epsilon;
    Verbose = options_.csminwel.verbosity;
    Save_files = options_.csminwel.Save_files;
    % Change some options.
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                nit = options_list{i,2};
              case 'InitialInverseHessian'
                H0 = eval(options_list{i,2});
              case 'TolFun'
                crit = options_list{i,2};
              case 'NumgradAlgorithm'
                numgrad = options_list{i,2};
              case 'NumgradEpsilon'
                epsilon = options_list{i,2};
              case 'verbosity'
                Verbose = options_list{i,2};
              case 'SaveFiles'
                Save_files = options_list{i,2};                
              otherwise
                warning(['csminwel: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    if options_.silent_optimizer
        Save_files = 0; 
        Verbose = 0;
    end    
    % Set flag for analytical gradient.
    if options_.analytic_derivation
        analytic_grad=1;
    else
        analytic_grad=[];
    end
    % Call csminwell.
    [fval,opt_par_values,grad,inverse_hessian_mat,itct,fcount,exitflag] = ...
        csminwel1(objective_function, start_par_value, H0, analytic_grad, crit, nit, numgrad, epsilon, Verbose, Save_files, varargin{:});
    hessian_mat=inv(inverse_hessian_mat);
  case 5
    if options_.analytic_derivation==-1 %set outside as code for use of analytic derivation
        analytic_grad=1;
        crit = options_.newrat.tolerance.f_analytic;
        newratflag = 0; %analytical Hessian
    else
        analytic_grad=0;
        crit=options_.newrat.tolerance.f;
        newratflag = options_.newrat.hess; %default
    end
    nit=options_.newrat.maxiter;
    Verbose = options_.newrat.verbosity;
    Save_files = options_.newrat.Save_files;
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                nit = options_list{i,2};
              case 'Hessian'
                flag=options_list{i,2};
                if options_.analytic_derivation && flag~=0
                    error('newrat: analytic_derivation is incompatible with numerical Hessian. Using analytic Hessian')
                else
                    newratflag=flag; 
                end
              case 'TolFun'
                crit = options_list{i,2};
              case 'verbosity'
                Verbose = options_list{i,2};
              case 'SaveFiles'
                Save_files = options_list{i,2};                
              otherwise
                warning(['newrat: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    if options_.silent_optimizer
        Save_files = 0; 
        Verbose = 0;
    end    
    hess_info.gstep=options_.gstep;
    hess_info.htol = 1.e-4;
    hess_info.h1=options_.gradient_epsilon*ones(n_params,1);
    [opt_par_values,hessian_mat,gg,fval,invhess,new_rat_hess_info] = newrat(objective_function,start_par_value,bounds,analytic_grad,crit,nit,0,Verbose, Save_files,hess_info,varargin{:});
    %hessian_mat is the plain outer product gradient Hessian
  case 6
    [opt_par_values, hessian_mat, Scale, fval] = gmhmaxlik(objective_function, start_par_value, ...
                                                      Initial_Hessian, options_.mh_jscale, bounds, prior_information.p2, options_.gmhmaxlik, options_.optim_opt, varargin{:});
  case 7
    % Matlab's simplex (Optimization toolbox needed).
    if isoctave && ~user_has_octave_forge_package('optim')
        error('Option mode_compute=7 requires the optim package')
    elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
        error('Option mode_compute=7 requires the Optimization Toolbox')
    end
    optim_options = optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6);
    if ~isempty(options_.optim_opt)
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if options_.silent_optimizer
        optim_options = optimset(optim_options,'display','off');
    end
    if ~isoctave
        [opt_par_values,fval,exitflag] = fminsearch(objective_function,start_par_value,optim_options,varargin{:});
    else
        % Under Octave, use a wrapper, since fminsearch() does not have a
        % 4th arg, and only has two output args
        func = @(x) objective_function(x,varargin{:});
        [opt_par_values,fval] = fminsearch(func,start_par_value,optim_options);
    end
  case 8
    % Dynare implementation of the simplex algorithm.
    simplexOptions = options_.simplex;
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                simplexOptions.maxiter = options_list{i,2};
              case 'TolFun'
                simplexOptions.tolerance.f = options_list{i,2};
              case 'TolX'
                simplexOptions.tolerance.x = options_list{i,2};
              case 'MaxFunEvals'
                simplexOptions.maxfcall = options_list{i,2};
              case 'MaxFunEvalFactor'
                simplexOptions.maxfcallfactor = options_list{i,2};
              case 'InitialSimplexSize'
                simplexOptions.delta_factor = options_list{i,2};
              case 'verbosity'
                simplexOptions.verbose = options_list{i,2};
              otherwise
                warning(['simplex: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    if options_.silent_optimizer
        simplexOptions.verbose = options_list{i,2};
    end    
    [opt_par_values,fval,exitflag] = simplex_optimization_routine(objective_function,start_par_value,simplexOptions,parameter_names,varargin{:});
  case 9
    % Set defaults
    H0 = (bounds(:,2)-bounds(:,1))*0.2;
    H0(~isfinite(H0)) = 0.01;
    cmaesOptions = options_.cmaes;
    cmaesOptions.LBounds = bounds(:,1);
    cmaesOptions.UBounds = bounds(:,2);
    % Modify defaults
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                cmaesOptions.MaxIter = options_list{i,2};
              case 'TolFun'
                cmaesOptions.TolFun = options_list{i,2};
              case 'TolX'
                cmaesOptions.TolX = options_list{i,2};
              case 'MaxFunEvals'
                cmaesOptions.MaxFunEvals = options_list{i,2};
              case 'verbosity'
                if options_list{i,2}==0
                    cmaesOptions.DispFinal  = 'off';   % display messages like initial and final message';
                    cmaesOptions.DispModulo = '0';   % [0:Inf], disp messages after every i-th iteration';
                end
              case 'SaveFiles'
                if options_list{i,2}==0
                  cmaesOptions.SaveVariables='off';
                  cmaesOptions.LogModulo = '0';    % [0:Inf] if >1 record data less frequently after gen=100';
                  cmaesOptions.LogTime   = '0';    % [0:100] max. percentage of time for recording data';
                end
              case 'CMAESResume'
                  if options_list{i,2}==1
                    cmaesOptions.Resume = 'yes';
                  end
              otherwise
                warning(['cmaes: Unknown option (' options_list{i,1}  ')!'])
            end
        end
    end
    if options_.silent_optimizer
        cmaesOptions.DispFinal  = 'off';   % display messages like initial and final message';
        cmaesOptions.DispModulo = '0';   % [0:Inf], disp messages after every i-th iteration';
        cmaesOptions.SaveVariables='off';
        cmaesOptions.LogModulo = '0';    % [0:Inf] if >1 record data less frequently after gen=100';
        cmaesOptions.LogTime   = '0';    % [0:100] max. percentage of time for recording data';
    end    
    warning('off','CMAES:NonfinitenessRange');
    warning('off','CMAES:InitialSigma');
    [x, fval, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes(func2str(objective_function),start_par_value,H0,cmaesOptions,varargin{:});
    opt_par_values=BESTEVER.x;
    fval=BESTEVER.f;
  case 10
    simpsaOptions = options_.simpsa;
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                simpsaOptions.MAX_ITER_TOTAL = options_list{i,2};
              case 'TolFun'
                simpsaOptions.TOLFUN = options_list{i,2};
              case 'TolX'
                tolx = options_list{i,2};
                if tolx<0
                    simpsaOptions = rmfield(simpsaOptions,'TOLX'); % Let simpsa choose the default.
                else
                    simpsaOptions.TOLX = tolx;
                end
              case 'EndTemparature'
                simpsaOptions.TEMP_END = options_list{i,2};
              case 'MaxFunEvals'
                simpsaOptions.MAX_FUN_EVALS = options_list{i,2};
              case 'verbosity'
                  if options_list{i,2} == 0
                    simpsaOptions.DISPLAY = 'none';
                  else
                    simpsaOptions.DISPLAY = 'iter';
                  end                      
              otherwise
                warning(['simpsa: Unknown option (' options_list{i,1}  ')!'])
            end
        end
    end
    if options_.silent_optimizer
        simpsaOptions.DISPLAY = 'none';
    end
    simpsaOptionsList = options2cell(simpsaOptions);
    simpsaOptions = simpsaset(simpsaOptionsList{:});
    [LB, UB]=set_bounds_to_finite_values(bounds, options_.huge_number);
    [opt_par_values, fval, exitflag] = simpsa(func2str(objective_function),start_par_value,LB,UB,simpsaOptions,varargin{:});
  case 11
    options_.cova_compute = 0;
    [opt_par_values, stdh, lb_95, ub_95, med_param] = online_auxiliary_filter(start_par_value, varargin{:});
  case 12
    [LB, UB] = set_bounds_to_finite_values(bounds, options_.huge_number);
    tmp = transpose([fieldnames(options_.particleswarm), struct2cell(options_.particleswarm)]);
    particleswarmOptions = optimoptions(@particleswarm);
    particleswarmOptions = optimoptions(particleswarmOptions, tmp{:});
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        SupportedListOfOptions = {'CreationFcn', 'Display', 'DisplayInterval', 'FunctionTolerance', ...
                                    'FunValCheck', 'HybridFcn', 'InertiaRange', 'InitialSwarmMatrix', 'InitialSwarmSpan', ...
                                    'MaxIterations', 'MaxStallIterations', 'MaxStallTime', 'MaxTime', ...
                                    'MinNeighborsFraction', 'ObjectiveLimit', 'OutputFcn', 'PlotFcn', 'SelfAdjustmentWeight', ...
                                    'SocialAdjustmentWeight', 'SwarmSize', 'UseParallel', 'UseVectorized'};
        for i=1:rows(options_list)
            if ismember(options_list{i,1}, SupportedListOfOptions)
                particleswarmOptions = optimoptions(particleswarmOptions, options_list{i,1}, options_list{i,2});
            else
                warning(['particleswarm: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    % Get number of instruments.
    numberofvariables = length(start_par_value);
    % Set objective function.
    objfun = @(x) objective_function(x, varargin{:});
    if ischar(particleswarmOptions.SwarmSize)
        eval(['particleswarmOptions.SwarmSize = ' particleswarmOptions.SwarmSize ';'])
    end
    if isempty(particleswarmOptions.InitialSwarmMatrix)
        particleswarmOptions.InitialSwarmMatrix = zeros(particleswarmOptions.SwarmSize, numberofvariables);
        p = 1;
        FVALS = zeros(particleswarmOptions.SwarmSize, 1);
        while p<=particleswarmOptions.SwarmSize
            candidate = rand(numberofvariables, 1).*(UB-LB)+LB;
            [fval, info, exit_flag] = objfun(candidate);
            if exit_flag
                particleswarmOptions.InitialSwarmMatrix(p,:) = transpose(candidate);
                FVALS(p) = fval;
                p = p + 1;
            end
        end
    end
    % Set penalty to the worst value of the objective function.
    TMP = [particleswarmOptions.InitialSwarmMatrix, FVALS];
    TMP = sortrows(TMP, length(start_par_value)+1);
    penalty = TMP(end,end);
    % Define penalized objective.
    objfun = @(x) penalty_objective_function(x, objective_function, penalty, varargin{:});
    % Minimize the penalized objective (note that the penalty is not updated).
    [opt_par_values, fval, exitflag, output] = particleswarm(objfun, length(start_par_value), LB, UB, particleswarmOptions);
    opt_par_values = opt_par_values(:);
  case 101
    solveoptoptions = options_.solveopt;
    if ~isempty(options_.optim_opt)
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'TolX'
                solveoptoptions.TolX = options_list{i,2};
              case 'TolFun'
                solveoptoptions.TolFun = options_list{i,2};
              case 'MaxIter'
                solveoptoptions.MaxIter = options_list{i,2};
              case 'verbosity'
                solveoptoptions.verbosity = options_list{i,2};
              case 'SpaceDilation'
                solveoptoptions.SpaceDilation = options_list{i,2};
              case 'LBGradientStep'
                solveoptoptions.LBGradientStep = options_list{i,2};
              otherwise
                warning(['solveopt: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    if options_.silent_optimizer
        solveoptoptions.verbosity = 0;
    end
    [opt_par_values,fval]=solvopt(start_par_value,objective_function,[],[],[],solveoptoptions,varargin{:});
  case 102
    if isoctave
        error('Optimization algorithm 2 is not available under Octave')
    elseif ~user_has_matlab_license('GADS_Toolbox')
        error('Optimization algorithm 2 requires the Global Optimization Toolbox')
    end
    % Set default optimization options for simulannealbnd.
    optim_options = saoptimset('display','iter','TolFun',1e-8);
    if ~isempty(options_.optim_opt)
        eval(['optim_options = saoptimset(optim_options,' options_.optim_opt ');']);
    end
    if options_.silent_optimizer
        optim_options = optimset(optim_options,'display','off');
    end
    func = @(x)objective_function(x,varargin{:});
    [opt_par_values,fval,exitflag,output] = simulannealbnd(func,start_par_value,bounds(:,1),bounds(:,2),optim_options);
  otherwise
    if ischar(minimizer_algorithm)
        if exist(minimizer_algorithm)
            [opt_par_values, fval, exitflag] = feval(minimizer_algorithm,objective_function,start_par_value,varargin{:});
        else
            error('No minimizer with the provided name detected.')
        end
    else
        error(['Optimization algorithm ' int2str(minimizer_algorithm) ' is unknown!'])
    end
end

end

function [LB, UB]=set_bounds_to_finite_values(bounds, huge_number)
    LB=bounds(:,1);
    LB(isinf(LB))=-huge_number;
    UB=bounds(:,2);
    UB(isinf(UB))=huge_number;
end
