function generate_trace_plots(chain_number)
%function generate_trace_plots(chain_number)
% Generates trace plots for all estimated parameters and the posterior
%
% INPUTS
%    chain_number:  [scalar]    number of the chain for which to construct the trace plots
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2016 Dynare Team
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

global M_ options_ estim_params_


% Get informations about the posterior draws:
MetropolisFolder = CheckPath('metropolis',M_.dname);
load_last_mh_history_file(MetropolisFolder, M_.fname);
if chain_number>record.Nblck
    error('generate_trace_plots:: chain number is bigger than existing number of chains')
end

trace_plot(options_,M_,estim_params_,'PosteriorDensity',chain_number)

for ii=1:size(estim_params_.param_vals,1)
    parameter_name=deblank(M_.param_names(estim_params_.param_vals(ii,1),:));
    trace_plot(options_,M_,estim_params_,'DeepParameter',chain_number,parameter_name)
end

for ii=1:size(estim_params_.var_exo,1)
    parameter_name=deblank(M_.exo_names(estim_params_.var_exo(ii,1),:));
    trace_plot(options_,M_,estim_params_,'StructuralShock',chain_number,parameter_name)
end

for ii=1:size(estim_params_.corrn,1)
    parameter_name_1=deblank(M_.endo_names(estim_params_.corrn(ii,1),:));
    parameter_name_2=deblank(M_.endo_names(estim_params_.corrn(ii,2),:));
    trace_plot(options_,M_,estim_params_,'MeasurementError',chain_number,parameter_name_1,parameter_name_2)
end

for ii=1:size(estim_params_.corrx,1)
    parameter_name_1=deblank(M_.exo_names(estim_params_.corrx(ii,1),:));
    parameter_name_2=deblank(M_.exo_names(estim_params_.corrx(ii,2),:));
    trace_plot(options_,M_,estim_params_,'StructuralShock',chain_number,parameter_name_1,parameter_name_2)
end