function  [par, logpost, accepted, neval] = posterior_sampler_iteration(TargetFun,last_draw, last_posterior, sampler_options,varargin)

% function [par, logpost, accepted, neval] = posterior_sampler_iteration(TargetFun,last_draw, last_posterior, sampler_options,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_)
% posterior samplers
%
% INPUTS
%   posterior_sampler_options:       posterior sampler options
%   options_:       structure storing the options

% OUTPUTS
%   posterior_sampler_options:       checked posterior sampler options
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2015 Dynare Team
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


posterior_sampling_method = sampler_options.posterior_sampling_method;
mh_bounds = sampler_options.bounds;

switch posterior_sampling_method
    case 'slice'
        
        [par, logpost, neval] = slice_sampler(TargetFun,last_draw, [mh_bounds.lb mh_bounds.ub], sampler_options,varargin{:});
        accepted = 1;
    case 'random_walk_metropolis_hastings'
        neval = 1;
        ProposalFun = sampler_options.proposal_distribution;
        proposal_covariance_Cholesky_decomposition = sampler_options.proposal_covariance_Cholesky_decomposition;
        n = sampler_options.n;

        par = feval(ProposalFun, last_draw, proposal_covariance_Cholesky_decomposition, n);
        if all( par(:) > mh_bounds.lb ) && all( par(:) < mh_bounds.ub )
            try
                logpost = - feval(TargetFun, par(:),varargin{:});
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
        r = logpost-last_posterior;
        if (logpost > -inf) && (log(rand) < r)
            accepted = 1;
        else
            accepted = 0;
            par = last_draw;
            logpost = last_posterior;
        end
    case 'independent_metropolis_hastings'
        neval = 1;
        ProposalFun = sampler_options.proposal_distribution;
        ProposalDensity = sampler_options.ProposalDensity;
        proposal_covariance_Cholesky_decomposition = sampler_options.proposal_covariance_Cholesky_decomposition;
        n = sampler_options.n;
        xparam1 = sampler_options.xparam1;
        par = feval(ProposalFun, sampler_options.xparam1, proposal_covariance_Cholesky_decomposition, n);
        if all( par(:) > mh_bounds.lb ) && all( par(:) < mh_bounds.ub )
            try
                logpost = - feval(TargetFun, par(:),varargin{:});
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
            r = logpost - last_posterior + ...
                log(feval(ProposalDensity, last_draw, xparam1, proposal_covariance_Cholesky_decomposition, n)) - ...
                log(feval(ProposalDensity, par, xparam1, proposal_covariance_Cholesky_decomposition, n));
        if (logpost > -inf) && (log(rand) < r)
            accepted = 1;
        else
            accepted = 0;
            par = last_draw;
            logpost = last_posterior;
        end
end