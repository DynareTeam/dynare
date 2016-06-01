function [fval,exit_flag,arg1,arg2] = penalty_objective_function(x0,fcn,penalty,varargin)
    [fval,info,exit_flag,arg1,arg2] = fcn(x0,varargin{:});
    
    if info(1) ~= 0
        fval = penalty + info(4);
    end
end