function [LIK,lik] = gaussian_filter(ReducedForm,Y,start,DynareOptions)
% Evaluates the likelihood of a non-linear model approximating the
% predictive (prior) and filtered (posterior) densities for state variables
% by gaussian distributions.
% Gaussian approximation is done by:
% - a Kronrod-Paterson gaussian quadrature with a limited number of nodes.
% Multidimensional quadrature is obtained by the Smolyak operator (ref: Winschel & Kratzig, 2010).
% - a spherical-radial cubature (ref: Arasaratnam & Haykin, 2008,2009).
% - a scaled unscented transform cubature (ref: )
% - Monte-Carlo draws from a multivariate gaussian distribution.
% First and second moments of prior and posterior state densities are computed
% from the resulting nodes/particles and allows to generate new distributions at the
% following observation.
% => The use of nodes is much faster than Monte-Carlo Gaussian particle and standard particles
% filters since it treats a lesser number of particles and there is no need
% of resampling.
% However, estimations may reveal biaised if the model is truly non-gaussian
% since predictive and filtered densities are unimodal.
%
% INPUTS
%    reduced_form_model     [structure] Matlab's structure describing the reduced form model.
%                                       reduced_form_model.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       reduced_form_model.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       reduced_form_model.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    smolyak_accuracy       [integer]   scalar.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.
% Copyright (C) 2009-2013 Dynare Team
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

persistent init_flag mf0 mf1
persistent nodes2 weights2 weights_c2 number_of_particles
persistent sample_size number_of_state_variables number_of_observed_variables

% Set default
if isempty(start)
    start = 1;
end

% Set persistent variables.
if isempty(init_flag)
    mf0 = ReducedForm.mf0;
    mf1 = ReducedForm.mf1;
    sample_size = size(Y,2);
    number_of_state_variables = length(mf0);
    number_of_observed_variables = length(mf1);
    number_of_particles = DynareOptions.particle.number_of_particles;
    init_flag = 1;
end

% compute gaussian quadrature nodes and weights on states and shocks
if isempty(nodes2)
    if DynareOptions.particle.distribution_approximation.cubature
        [nodes2,weights2] = spherical_radial_sigma_points(number_of_state_variables);
        weights_c2 = weights2;
    elseif DynareOptions.particle.distribution_approximation.unscented
        [nodes2,weights2,weights_c2] = unscented_sigma_points(number_of_state_variables,DynareOptions);
    else
        if ~DynareOptions.particle.distribution_approximation.montecarlo
            error('Estimation: This approximation for the proposal is not implemented or unknown!')
        end
    end
end

if DynareOptions.particle.distribution_approximation.montecarlo
    set_dynare_seed('default');
end

% Get covariance matrices
Q = ReducedForm.Q;
H = ReducedForm.H;
if isempty(H)
    H = 0;
    H_lower_triangular_cholesky = 0;
else
    H_lower_triangular_cholesky = reduced_rank_cholesky(H)';
end

% Get initial condition for the state vector.
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);
Q_lower_triangular_cholesky = reduced_rank_cholesky(Q)';

% Initialization of the likelihood.
const_lik = (2*pi)^(number_of_observed_variables/2) ;
lik  = NaN(sample_size,1);
LIK  = NaN;

SampleWeights = 1/number_of_particles ;
ks = 0 ;
%Estimate = zeros(number_of_state_variables,sample_size) ;
%V_Estimate = zeros(number_of_state_variables,number_of_state_variables,sample_size) ;
for t=1:sample_size
    % build the proposal
    [PredictedStateMean,PredictedStateVarianceSquareRoot,StateVectorMean,StateVectorVarianceSquareRoot] = ...
        gaussian_filter_bank(ReducedForm,Y(:,t),StateVectorMean,StateVectorVarianceSquareRoot,Q_lower_triangular_cholesky,H_lower_triangular_cholesky,H,DynareOptions) ;
    %Estimate(:,t) = PredictedStateMean ;
    %V_Estimate(:,:,t) = PredictedStateVarianceSquareRoot ;
    if DynareOptions.particle.distribution_approximation.cubature || DynareOptions.particle.distribution_approximation.unscented
        StateParticles = bsxfun(@plus,StateVectorMean,StateVectorVarianceSquareRoot*nodes2') ;
        IncrementalWeights = ...
                    gaussian_densities(Y(:,t),StateVectorMean,...
                                        StateVectorVarianceSquareRoot,PredictedStateMean,...
                                        PredictedStateVarianceSquareRoot,StateParticles,H,const_lik,...
                                        weights2,weights_c2,ReducedForm,DynareOptions) ;
        SampleWeights = weights2.*IncrementalWeights ;
        SumSampleWeights = sum(SampleWeights) ;
        lik(t) = log(SumSampleWeights) ;
        SampleWeights = SampleWeights./SumSampleWeights ;
    else % Monte-Carlo draws
        StateParticles = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),StateVectorMean) ;
        IncrementalWeights = ...
                    gaussian_densities(Y(:,t),StateVectorMean,...
                                        StateVectorVarianceSquareRoot,PredictedStateMean,...
                                        PredictedStateVarianceSquareRoot,StateParticles,H,const_lik,...
                                        1/number_of_particles,1/number_of_particles,ReducedForm,DynareOptions) ;
        SampleWeights = SampleWeights.*IncrementalWeights ;
        SumSampleWeights = sum(SampleWeights) ;
        %VarSampleWeights = IncrementalWeights-SumSampleWeights ;
        %VarSampleWeights = VarSampleWeights*VarSampleWeights'/(number_of_particles-1) ;
        lik(t) = log(SumSampleWeights) ; %+ .5*VarSampleWeights/(number_of_particles*(SumSampleWeights*SumSampleWeights)) ;
        SampleWeights = SampleWeights./SumSampleWeights ;
        Neff = 1/sum(bsxfun(@power,SampleWeights,2)) ;
        if (Neff<.5*sample_size && DynareOptions.particle.resampling.status.generic) || DynareOptions.particle.resampling.status.systematic
            ks = ks + 1 ;
            StateParticles = resample(StateParticles',SampleWeights,DynareOptions)' ;
            StateVectorMean = mean(StateParticles,2) ;
            StateVectorVarianceSquareRoot = reduced_rank_cholesky( (StateParticles*StateParticles')/(number_of_particles-1) - StateVectorMean*(StateVectorMean') )';
            SampleWeights = 1/number_of_particles ;
        elseif DynareOptions.particle.resampling.status.none
            StateVectorMean = (sampleWeights*StateParticles)' ;
            temp = sqrt(SampleWeights').*StateParticles ;
            StateVectorVarianceSquareRoot = reduced_rank_cholesky( temp'*temp - StateVectorMean*(StateVectorMean') )';
        end
    end
end

LIK = -sum(lik(start:end));
