function [x,info,fvec,fjac] = dynare_solve(func,x,options,varargin)
% function [x,info,fvec,fjac] = dynare_solve(func,x,options,varargin)
% proposes different solvers
%
% INPUTS
%    func:             name of the function to be solved
%    x:                guess values
%    options:          struct of Dynare options
%    varargin:         list of arguments following jacobian_flag
%
% OUTPUTS
%    x:                solution
%    info=1:           the model can not be solved
%    fvec:             Function value (used for debugging when check=1)
%    fjac:             Jacobian (used for debugging when check=1)
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2016 Dynare Team
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

% jacobian_flag=1:  jacobian given by the 'func' function
% jacobian_flag=0:  jacobian obtained numerically
jacobian_flag = options.jacobian_flag;

% Set tolerance parameter depending the the caller function.
stack = dbstack;
if isoctave
    [pathstr,name,ext]=fileparts(stack(2).file);
    caller_file_name=[name,ext];
else
    caller_file_name=stack(2).file;
end
if strcmp(caller_file_name,'simulation_core.m') || strcmp(caller_file_name,'solve_stacked_problem.m')
    tolf = options.dynatol.f;
else
    tolf = options.solve_tolf;
end

if strcmp(caller_file_name,'dyn_ramsey_static.m')
    maxit = options.ramsey.maxit;
else
    maxit = options.steady.maxit;
end


info = 0;
nn = size(x,1);

% checking initial values
if jacobian_flag
    [fvec,fjac] = feval(func,x,varargin{:});
    if any(any(isinf(fjac) | isnan(fjac)))
        info=1;
        x = NaN(size(fvec));
        return
    end
else
    fvec = feval(func,x,varargin{:});
    fjac = zeros(nn,nn) ;
end

i = find(~isfinite(fvec));

if ~isempty(i)
    info = 1;
    x = NaN(size(fvec));
    return;
end

% this test doesn't check complementarity conditions and is not used for
% mixed complementarity problems
if (~ismember(options.solve_algo,[10,11])) && (max(abs(fvec)) < tolf)
    return ;
end

if options.solve_algo == 0
    if ~isoctave
        if ~user_has_matlab_license('optimization_toolbox')
            error('You can''t use solve_algo=0 since you don''t have MATLAB''s Optimization Toolbox')
        end
    end
    options4fsolve=optimset('fsolve');
    options4fsolve.MaxFunEvals = 50000;
    options4fsolve.MaxIter = maxit;
    options4fsolve.TolFun = tolf;
    if options.debug==1
        options4fsolve.Display = 'final';
    else
        options4fsolve.Display = 'off';
    end
    if jacobian_flag
        options4fsolve.Jacobian = 'on';
    else
        options4fsolve.Jacobian = 'off';
    end
    if ~isoctave
        [x,fval,exitval,output] = fsolve(func,x,options4fsolve,varargin{:});
    else
        % Under Octave, use a wrapper, since fsolve() does not have a 4th arg
        if ischar(func)
            func2 = str2func(func);
        else
            func2 = func;
        end
        func = @(x) func2(x, varargin{:});
        % The Octave version of fsolve does not converge when it starts from the solution
        fvec = feval(func,x);
        if max(abs(fvec)) >= tolf
            [x,fval,exitval,output] = fsolve(func,x,options4fsolve);
        else
            exitval = 3;
        end;
    end

    if exitval == 1
        info = 0;
    elseif exitval > 1
        if ischar(func)
            func2 = str2func(func);
        else
            func2 = func;
        end
        func = @(x) func2(x, varargin{:});
        fvec = feval(func,x);
        if max(abs(fvec)) >= tolf
            info = 1;
        else
            info = 0;        
        end
    else
        info = 1;
    end
elseif options.solve_algo == 1
        [x,info]=solve1(func,x,1:nn,1:nn,jacobian_flag,options.gstep, ...
                    tolf,options.solve_tolx, ...
                    maxit,options.debug,varargin{:});
elseif options.solve_algo == 9
        [x,info]=trust_region(func,x,1:nn,1:nn,jacobian_flag,options.gstep, ...
                    tolf,options.solve_tolx, ...
                    maxit,options.debug,varargin{:});
elseif options.solve_algo == 2 || options.solve_algo == 4

    if options.solve_algo == 2
        solver = @solve1;
    else
        solver = @trust_region;
    end

    if ~jacobian_flag
        fjac = zeros(nn,nn) ;
        dh = max(abs(x),options.gstep(1)*ones(nn,1))*eps^(1/3);
        for j = 1:nn
            xdh = x ;
            xdh(j) = xdh(j)+dh(j) ;
            fjac(:,j) = (feval(func,xdh,varargin{:}) - fvec)./dh(j) ;
        end
    end

    [j1,j2,r,s] = dmperm(fjac);

    if options.debug
        disp(['DYNARE_SOLVE (solve_algo=2|4): number of blocks = ' num2str(length(r))]);
    end

    for i=length(r)-1:-1:1
        if options.debug
            disp(['DYNARE_SOLVE (solve_algo=2|4): solving block ' num2str(i) ', of size ' num2str(r(i+1)-r(i)) ]);
        end
        [x,info]=solver(func,x,j1(r(i):r(i+1)-1),j2(r(i):r(i+1)-1),jacobian_flag, ...
                        options.gstep, ...
                        tolf,options.solve_tolx, ...
                        maxit,options.debug,varargin{:});
        if info
            return
        end
    end
    fvec = feval(func,x,varargin{:});
    if max(abs(fvec)) > tolf
        [x,info]=solver(func,x,1:nn,1:nn,jacobian_flag, ...
                        options.gstep, tolf,options.solve_tolx, ...
                        maxit,options.debug,varargin{:});
    end
elseif options.solve_algo == 3
    if jacobian_flag
        [x,info] = csolve(func,x,func,1e-6,500,varargin{:});
    else
        [x,info] = csolve(func,x,[],1e-6,500,varargin{:});
    end
elseif options.solve_algo == 10
    % LMMCP
    olmmcp = options.lmmcp;

    [x,fval,exitflag] = lmmcp(func,x,olmmcp.lb,olmmcp.ub,olmmcp,varargin{:});
    if exitflag == 1
        info = 0;
    else
        info = 1;
    end
elseif options.solve_algo == 11
    % PATH mixed complementary problem
    % PATH linear mixed complementary problem
    if ~exist('mcppath')
        error(['PATH can''t be provided with Dynare. You need to install it ' ...
               'yourself and add its location to Matlab/Octave path before ' ...
               'running Dynare'])
    end
    omcppath = options.mcppath;
    global mcp_data
    mcp_data.func = func;
    mcp_data.args = varargin;
    info=0;
    try
        [x,fval,jac,mu] = pathmcp(x,omcppath.lb,omcppath.ub,'mcp_func',omcppath.A,omcppath.b,omcppath.t,omcppath.mu0);
    catch
        info = 1;
    end
else
    error('DYNARE_SOLVE: option solve_algo must be one of [0,1,2,3,4,9,10,11]')
end

