function [ProposalStateVector,Weights] = conditional_filter_proposal(ReducedForm,obs,StateVectors,SampleWeights,Q_lower_triangular_cholesky,H_lower_triangular_cholesky,H,DynareOptions,normconst2) 
%
% Computes the proposal for each past particle using Gaussian approximations 
% for the state errors and the Kalman filter 
% 
% INPUTS
%    reduced_form_model     [structure] Matlab's structure describing the reduced form model.
%                                       reduced_form_model.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       reduced_form_model.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       reduced_form_model.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.
% Copyright (C) 2012-2013 Dynare Team
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
%
% AUTHOR(S) frederic DOT karame AT univ DASH lemans DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

persistent init_flag2 mf0 mf1
persistent number_of_state_variables number_of_observed_variables 
persistent number_of_structural_innovations 

% Set local state space model (first-order approximation).
ghx  = ReducedForm.ghx;
ghu  = ReducedForm.ghu;
% Set local state space model (second-order approximation).
ghxx = ReducedForm.ghxx;
ghuu = ReducedForm.ghuu;
ghxu = ReducedForm.ghxu;

if any(any(isnan(ghx))) || any(any(isnan(ghu))) || any(any(isnan(ghxx))) || any(any(isnan(ghuu))) || any(any(isnan(ghxu))) || ...
    any(any(isinf(ghx))) || any(any(isinf(ghu))) || any(any(isinf(ghxx))) || any(any(isinf(ghuu))) || any(any(isinf(ghxu))) ...
    any(any(abs(ghx)>1e4)) || any(any(abs(ghu)>1e4)) || any(any(abs(ghxx)>1e4)) || any(any(abs(ghuu)>1e4)) || any(any(abs(ghxu)>1e4))
    ghx
    ghu
    ghxx
    ghuu
    ghxu
end

constant = ReducedForm.constant;
state_variables_steady_state = ReducedForm.state_variables_steady_state;

% Set persistent variables.
if isempty(init_flag2)
    mf0 = ReducedForm.mf0;
    mf1 = ReducedForm.mf1;
    number_of_state_variables = length(mf0);
    number_of_observed_variables = length(mf1);
    number_of_structural_innovations = length(ReducedForm.Q);
    init_flag2 = 1;
end

if strcmpi(DynareOptions.particle.IS_approximation_method,'cubature') || strcmpi(DynareOptions.particle.IS_approximation_method,'monte-carlo')
    [nodes,weights] = spherical_radial_sigma_points(number_of_structural_innovations) ; 
    weights_c = weights ;
end    
if strcmpi(DynareOptions.particle.IS_approximation_method,'quadrature')
    [nodes,weights] = nwspgr('GQN',number_of_structural_innovations,DynareOptions.particle.smolyak_accuracy) ;
    weights_c = weights ;
end    
if strcmpi(DynareOptions.particle.IS_approximation_method,'unscented')
    [nodes,weights,weights_c] = unscented_sigma_points(number_of_structural_innovations,DynareOptions) ; 
end

epsilon = Q_lower_triangular_cholesky*(nodes') ;
yhat = repmat(StateVectors-state_variables_steady_state,1,size(epsilon,2)) ;

tmp = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,DynareOptions.threads.local_state_space_iteration_2);

PredictedStateMean = tmp(mf0,:)*weights ;
PredictedObservedMean = tmp(mf1,:)*weights;

if strcmpi(DynareOptions.particle.IS_approximation_method,'cubature') || ...
   strcmpi(DynareOptions.particle.IS_approximation_method,'monte-carlo')
    PredictedStateMean = sum(PredictedStateMean,2) ;
    PredictedObservedMean = sum(PredictedObservedMean,2) ;
    dState = bsxfun(@minus,tmp(mf0,:),PredictedStateMean)'.*sqrt(weights) ;
    dObserved = bsxfun(@minus,tmp(mf1,:),PredictedObservedMean)'.*sqrt(weights);
    big_mat = [dObserved  dState ; [H_lower_triangular_cholesky zeros(number_of_observed_variables,number_of_state_variables)] ] ;
    [mat1,mat] = qr2(big_mat,0) ;
    mat = mat' ; 
    clear('mat1');
    PredictedObservedVarianceSquareRoot = mat(1:number_of_observed_variables,1:number_of_observed_variables) ;
    CovarianceObservedStateSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),1:number_of_observed_variables) ;
    StateVectorVarianceSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),number_of_observed_variables+(1:number_of_state_variables)) ;
    StateVectorMean = PredictedStateMean + (CovarianceObservedStateSquareRoot/PredictedObservedVarianceSquareRoot)*(obs - PredictedObservedMean) ; 
end 
if strcmpi(DynareOptions.particle.IS_approximation_method,'quadrature') || ... 
   strcmpi(DynareOptions.particle.IS_approximation_method,'unscented') 
    dState = bsxfun(@minus,tmp(mf0,:),PredictedStateMean);
    dObserved = bsxfun(@minus,tmp(mf1,:),PredictedObservedMean);
    PredictedStateVariance = dState*diag(weights_c)*dState';
    PredictedObservedVariance = dObserved*diag(weights_c)*dObserved' + H;
    PredictedStateAndObservedCovariance = dState*diag(weights_c)*dObserved';
    KalmanFilterGain = PredictedStateAndObservedCovariance/PredictedObservedVariance ; 
    StateVectorMean = PredictedStateMean + KalmanFilterGain*(obs - PredictedObservedMean);
    StateVectorVariance = PredictedStateVariance - KalmanFilterGain*PredictedObservedVariance*KalmanFilterGain';
    StateVectorVariance = .5*(StateVectorVariance+StateVectorVariance');
    StateVectorVarianceSquareRoot = reduced_rank_cholesky(StateVectorVariance)'; 
end 
    
ProposalStateVector = StateVectorVarianceSquareRoot*randn(size(StateVectorVarianceSquareRoot,2),1)+StateVectorMean ;
ypred = measurement_equations(ProposalStateVector,ReducedForm,DynareOptions) ;
foo = H_lower_triangular_cholesky \ (obs - ypred) ;
likelihood = exp(-0.5*(foo)'*foo)/normconst2 + 1e-99 ;			
Weights = SampleWeights.*likelihood ;	

