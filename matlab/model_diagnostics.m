function model_diagnostics(M,options,oo)
% function model_diagnostics(M,options,oo)
%   computes various diagnostics on the model 
% INPUTS
%   M         [matlab structure] Definition of the model.           
%   options   [matlab structure] Global options.
%   oo        [matlab structure] Results 
%    
% OUTPUTS
%   none
%    
% ALGORITHM
%   ...
%    
% SPECIAL REQUIREMENTS
%   none.
%  

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

global jacob

endo_nbr = M.endo_nbr;
endo_names = M.endo_names;
lead_lag_incidence = M.lead_lag_incidence;
maximum_endo_lag = M.maximum_endo_lag;

%
% missing variables at the current period
%
k = find(lead_lag_incidence(maximum_endo_lag+1,:)==0);
if ~isempty(k)
    disp(['The following endogenous variables aren''t present at ' ...
          'the current period in the model:'])
    for i=1:length(k)
        disp(endo_names(k(i),:))
    end
end

%
% check steady state
%
info = 0;

if M.exo_nbr == 0
    oo.exo_steady_state = [] ;
end

% check if ys is steady state
[dr.ys,params,check1]=evaluate_steady_state(oo.steady_state,M,options,oo,1);

% testing for problem
if check1(1)
    disp('model diagnostic can''t obtain the steady state')
    if any(isnan(dr.ys))
        disp(['model diagnostic obtains a steady state with NaNs'])
    end
    if any(isinf(dr.ys))
        disp(['model diagnostic obtains a steady state with Inf'])
    end
    return;
end

if ~isreal(dr.ys)
    disp(['model diagnostic obtains a steady state with complex ' ...
          'numbers'])
    return
end

%
% singular Jacobian of static model
%
if ~isfield(M,'block_structure_stat')
    nb = 1;
else
    nb = length(M.block_structure_stat.block);
end

exo = [oo.exo_steady_state; oo.exo_det_steady_state];
for b=1:nb
    if options.bytecode
        if nb == 1
            [chk, res, jacob] = bytecode(dr.ys, exo, M.params, dr.ys, 1, exo, ...
                                    'evaluate', 'static');
        else
            [chk, res, jacob] = bytecode(dr.ys, exo, M.params, dr.ys, 1, exo, ...
                                         'evaluate', 'static',['block=' ...
                                int2str(b)]);
        end
    else
        [res,jacob]=feval([M.fname '_static'],dr.ys,exo,M.params);
    end
    rank_jacob = rank(jacob);
    if rank_jacob < size(jacob,1)
        disp(['model_diagnostic: the Jacobian of the static model is ' ...
              'singular'])
        disp(['there is ' num2str(endo_nbr-rank_jacob) ...
              ' colinear relationships between the variables and the equations'])
        ncol = null(jacob);
        n_rel = size(ncol,2);
        for i = 1:n_rel
            if n_rel  > 1
                disp(['Relation ' int2str(i)])
            end
            disp('Colinear variables:')
            for j=1:10
                k = find(abs(ncol(:,i)) > 10^-j);
                if max(abs(jacob(:,k)*ncol(k,i))) < 1e-6
                    break
                end
            end
            disp(endo_names(k,:))
        end
        neq = null(jacob');
        n_rel = size(neq,2);
        for i = 1:n_rel
            if n_rel  > 1
                disp(['Relation ' int2str(i)])
            end
            disp('Colinear equations')
            for j=1:10
                k = find(abs(neq(:,i)) > 10^-j);
                if max(abs(jacob(k,:)'*neq(k,i))) < 1e-6
                    break
                end
            end
            disp(k')
        end
    end
end

