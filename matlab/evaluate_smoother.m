function oo = evaluate_smoother(parameters)
% Evaluate the smoother at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean') or a vector of values for
%                  the (estimated) parameters of the model.
%
%
% OUTPUTS
%    o oo       [structure]  results:
%                              - SmoothedVariables
%                              - SmoothedShocks
%                              - SmoothedVariables
%                              - SmoothedVariables
%                              - SmoothedVariables
%                              - SmoothedVariables
%                              - SmoothedVariables
%                              - SmoothedVariables
%
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function use persistent variables for the dataset and the description of the missing observations. Consequently, if this function
%     is called more than once (by changing the value of parameters) the sample *must not* change.

% Copyright (C) 2010-2012 Dynare Team
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

global options_ M_ bayestopt_ oo_

persistent dataset_

if nargin==0
    parameters = 'posterior_mode';
end

if ischar(parameters)
    switch parameters
      case 'posterior_mode'
        parameters = get_posterior_parameters('mode');
      case 'posterior_mean'
        parameters = get_posterior_parameters('mean');
      case 'posterior_median'
        parameters = get_posterior_parameters('median');
      case 'prior_mode'
        parameters = bayestopt_.p5(:);
      case 'prior_mean'
        parameters = bayestopt_.p1;
      case 'calibration'
        if isempty(oo_.dr)
            error('You must run ''stoch_simul'' first.');
        end
        parameters = [];
      otherwise
        disp('evaluate_smoother:: If the input argument is a string, then it has to be equal to:')
        disp('                     ''posterior_mode'', ')
        disp('                     ''posterior_mean'', ')
        disp('                     ''posterior_median'', ')
        disp('                     ''prior_mode'' or')
        disp('                     ''prior_mean''.')
        disp('                     ''calibration''.')
        error
    end
end

if isempty(dataset_)
    % Load and transform data.
    transformation = [];
    if options_.loglinear && ~options_.logdata
        transformation = @log;
    end
    xls.sheet = options_.xls_sheet;
    xls.range = options_.xls_range;

    if ~isfield(options_,'nobs')
        options_.nobs = [];
    end

    dataset_ = initialize_dataset(options_.datafile,options_.varobs,options_.first_obs,options_.nobs,transformation,options_.prefilter,xls);
    options_.nobs = dataset_.info.ntobs;

    % Determine if a constant is needed.
    if options_.steadystate_flag% if the *_steadystate.m file is provided.
        [ys,params,info] = evaluate_steady_state(oo_.steady_state,M_,options_,oo_,1);
        if size(ys,1) < M_.endo_nbr
            if length(M_.aux_vars) > 0
                ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars,...
                                                            M_.fname,...
                                                            zeros(M_.exo_nbr,1),...
                                                            oo_.exo_det_steady_state,...
                                                            M_.params,...
                                                            options_.bytecode);
            else
                error([M_.fname '_steadystate.m doesn''t match the model']);
            end
        end
        oo_.steady_state = ys;
    else% if the steady state file is not provided.
        [dd,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
        oo_.steady_state = dd.ys; clear('dd');
    end
    if all(abs(oo_.steady_state(bayestopt_.mfys))<1e-9)
        options_.noconstant = 1;
    else
        options_.noconstant = 0;
    end
end

pshape_original   = bayestopt_.pshape;
bayestopt_.pshape = Inf(size(bayestopt_.pshape));
clear('priordens')

[atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = ...
    DsgeSmoother(parameters,dataset_.info.ntobs,dataset_.data,dataset_.missing.aindex,dataset_.missing.state);

oo.Smoother.SteadyState = ys;
oo.Smoother.TrendCoeffs = trend_coeff;
if options_.filter_covariance
    oo.Smoother.Variance = P;
end
i_endo = bayestopt_.smoother_saved_var_list;
if options_.nk ~= 0
    oo.FilteredVariablesKStepAhead = ...
        aK(options_.filter_step_ahead,i_endo,:);
    if ~isempty(PK)
        oo.FilteredVariablesKStepAheadVariances = ...
            PK(options_.filter_step_ahead,i_endo,i_endo,:);
    end
    if ~isempty(decomp)
        oo.FilteredVariablesShockDecomposition = ...
            decomp(options_.filter_step_ahead,i_endo,:,:);
    end
end
dr = oo_.dr;
order_var = oo_.dr.order_var;
for i=bayestopt_.smoother_saved_var_list'
    i1 = order_var(bayestopt_.smoother_var_list(i));
    eval(['oo.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ' = atT(i,:)'';']);
    eval(['oo.FilteredVariables.' deblank(M_.endo_names(i1,:)) ' = squeeze(aK(1,i,:));']);
    eval(['oo.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ' = updated_variables(i,:)'';']);
end
for i=1:M_.exo_nbr
    eval(['oo.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
end

oo.dr = oo_.dr;

bayestopt_.pshape = pshape_original;