function [xparam,std_param,lb_95,ub_95,median_param] = online_auxiliary_filter(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)

% Carvalho & Lopes particle filter = auxiliary particle filter including Liu & West filter on parameters.
%
% INPUTS
%    ReducedForm     [structure] Matlab's structure describing the reduced form model.
%                                       ReducedForm.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       ReducedForm.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       ReducedForm.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                     [integer]   pp*1 vector of indices.
%    number_of_particles    [integer]   scalar.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2013 Dynare Team
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
persistent Y init_flag mf0 mf1 bounds number_of_particles number_of_parameters liu_west_delta liu_west_chol_sigma_bar
persistent start_param sample_size number_of_observed_variables number_of_structural_innovations

% Set seed for randn().
%set_dynare_seed('mt19937ar',1234) ;
set_dynare_seed('default') ;
pruning = DynareOptions.particle.pruning;
second_resample = 1 ;
variance_update = 1 ;

% initialization of state particles
[ys,trend_coeff,exit_flag,info,Model,DynareOptions,BayesInfo,DynareResults,ReducedForm] = ... 
            solve_model_for_online_filter(1,xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults) ;

% Set persistent variables.
if isempty(init_flag)
    mf0 = ReducedForm.mf0;
    mf1 = ReducedForm.mf1;
    number_of_particles = DynareOptions.particle.number_of_particles;
    number_of_parameters = size(xparam1,1) ;
    Y = DynareDataset.data ;
    sample_size = size(Y,2);
    number_of_observed_variables = length(mf1);
    number_of_structural_innovations = length(ReducedForm.Q);
    liu_west_delta = DynareOptions.particle.liu_west_delta ;
    liu_west_chol_sigma_bar = DynareOptions.particle.liu_west_chol_sigma_bar*eye(number_of_parameters) ;
    %start_param = xparam1 ; 
    % Conditions initiales 
    %liu_west_chol_sigma_bar = bsxfun(@times,eye(number_of_parameters),BayesInfo.p2) ; 
    start_param = BayesInfo.p1 ;
    bounds = [BayesInfo.lb BayesInfo.ub] ;
    init_flag = 1;
end

% Get initial conditions for the state particles 
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),StateVectorMean);
if pruning
    StateVectors_ = StateVectors;
end
        
% parameters for the Liu & West filter 
h_square = (3*liu_west_delta-1)/(2*liu_west_delta) ;
h_square = 1-h_square*h_square ;
small_a = sqrt(1-h_square) ;

% Initialization of parameter particles 
xparam = zeros(number_of_parameters,number_of_particles) ;
stderr = sqrt(bsxfun(@power,bounds(:,2)+bounds(:,1),2)/12)/1000 ;
i = 1 ;
while i<=number_of_particles
    candidate = start_param + 10*liu_west_chol_sigma_bar*randn(number_of_parameters,1) ;
    %candidate = start_param + bsxfun(@times,stderr,randn(number_of_parameters,1)) ;
    if all(candidate(:) >= bounds(:,1)) && all(candidate(:) <= bounds(:,2))
        xparam(:,i) = candidate(:) ;
        i = i+1 ;
    end
end

%xparam = bsxfun(@plus,bounds(:,1),bsxfun(@times,(bounds(:,2)-bounds(:,1)),rand(number_of_parameters,number_of_particles))) ;

% Initialization of the weights of particles.
weights = ones(1,number_of_particles)/number_of_particles ;

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables;
mean_xparam = zeros(number_of_parameters,sample_size) ;
median_xparam = zeros(number_of_parameters,sample_size) ;
std_xparam = zeros(number_of_parameters,sample_size) ;
lb95_xparam = zeros(number_of_parameters,sample_size) ;
ub95_xparam = zeros(number_of_parameters,sample_size) ;

%% The Online filter 
for t=1:sample_size
    disp(t) 
    % Moments of parameters particles distribution 
    m_bar = xparam*(weights') ;
    temp = bsxfun(@minus,xparam,m_bar) ;
    sigma_bar = (bsxfun(@times,weights,temp))*(temp') ;
    if variance_update==1
        chol_sigma_bar = chol(h_square*sigma_bar)' ;
    end
    % Prediction (without shocks)
    ObservedVariables = zeros(number_of_observed_variables,number_of_particles) ;
    for i=1:number_of_particles
        % model resolution 
        [ys,trend_coeff,exit_flag,info,Model,DynareOptions,BayesInfo,DynareResults,ReducedForm] = ... 
            solve_model_for_online_filter(t,xparam(:,i),DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults) ;
        steadystate = ReducedForm.steadystate;
        state_variables_steady_state = ReducedForm.state_variables_steady_state;
        % Set local state space model (second-order approximation).
        constant = ReducedForm.constant;
        ghx  = ReducedForm.ghx;
        ghu  = ReducedForm.ghu;
        ghxx = ReducedForm.ghxx;
        ghuu = ReducedForm.ghuu;
        ghxu = ReducedForm.ghxu;
        % particle likelihood contribution  
        yhat = bsxfun(@minus,StateVectors(:,i),state_variables_steady_state);
        if pruning
            yhat_ = bsxfun(@minus,StateVectors_(:,i),state_variables_steady_state);
            [tmp, tmp_] = local_state_space_iteration_2(yhat,zeros(number_of_structural_innovations,1),ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,steadystate,DynareOptions.threads.local_state_space_iteration_2);
        else
            tmp = local_state_space_iteration_2(yhat,zeros(number_of_structural_innovations,1),ghx,ghu,constant,ghxx,ghuu,ghxu,DynareOptions.threads.local_state_space_iteration_2);
        end
        ObservedVariables(:,i) = tmp(mf1,:) ;
    end
    PredictedObservedMean = sum(bsxfun(@times,weights,ObservedVariables),2) ;
    PredictionError = bsxfun(@minus,Y(:,t),ObservedVariables);
    dPredictedObservedMean = bsxfun(@minus,ObservedVariables,PredictedObservedMean);
    PredictedObservedVariance = bsxfun(@times,weights,dPredictedObservedMean)*dPredictedObservedMean' + ReducedForm.H ;
    wtilde = exp(-.5*(const_lik+log(det(PredictedObservedVariance))+sum(PredictionError.*(PredictedObservedVariance\PredictionError),1))) ;
    % unormalized weights and observation likelihood contribution 
    tau_tilde = weights.*wtilde ;
    sum_tau_tilde = sum(tau_tilde) ;
    % particles selection 
    tau_tilde = tau_tilde/sum_tau_tilde ;
    indx = index_resample(0,tau_tilde',DynareOptions);
    StateVectors = StateVectors(:,indx) ;
    if pruning
        StateVectors_ = StateVectors_(:,indx) ;
    end
    xparam = bsxfun(@plus,(1-small_a).*m_bar,small_a.*xparam) ;
    xparam = xparam(:,indx) ;
    wtilde = wtilde(indx) ;
    % draw in the new distributions 
    i = 1 ;
    while i<=number_of_particles
        candidate = xparam(:,i) + chol_sigma_bar*randn(number_of_parameters,1) ;
        if all(candidate >= bounds(:,1)) && all(candidate <= bounds(:,2))
            xparam(:,i) = candidate ;
            % model resolution for new parameters particles
            [ys,trend_coeff,exit_flag,info,Model,DynareOptions,BayesInfo,DynareResults,ReducedForm] = ...
                solve_model_for_online_filter(t,xparam(:,i),DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults) ;
            steadystate = ReducedForm.steadystate;
            state_variables_steady_state = ReducedForm.state_variables_steady_state;
            % Set local state space model (second order approximation).
            constant = ReducedForm.constant;
            ghx  = ReducedForm.ghx;
            ghu  = ReducedForm.ghu;
            ghxx = ReducedForm.ghxx;
            ghuu = ReducedForm.ghuu;
            ghxu = ReducedForm.ghxu;
            % Get covariance matrices and structural shocks 
            epsilon = chol(ReducedForm.Q)'*randn(number_of_structural_innovations,1) ;
            % compute particles likelihood contribution 
            yhat = bsxfun(@minus,StateVectors(:,i),state_variables_steady_state);
            if pruning
                yhat_ = bsxfun(@minus,StateVectors_(:,i),state_variables_steady_state);
                [tmp, tmp_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,steadystate,DynareOptions.threads.local_state_space_iteration_2);
                StateVectors_(:,i) = tmp_(mf0,:) ;
            else
                tmp = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,DynareOptions.threads.local_state_space_iteration_2);
            end
            StateVectors(:,i) = tmp(mf0,:) ;
            ObservedVariables(:,i) = tmp(mf1,:) ;
            i = i+1 ;
        end
    end
    PredictedObservedMean = sum(bsxfun(@times,weights,ObservedVariables),2) ;
    PredictionError = bsxfun(@minus,Y(:,t),ObservedVariables);
    dPredictedObservedMean = bsxfun(@minus,ObservedVariables,PredictedObservedMean);
    PredictedObservedVariance = bsxfun(@times,weights,dPredictedObservedMean)*dPredictedObservedMean' + ReducedForm.H ;
    lnw = exp(-.5*(const_lik+log(det(PredictedObservedVariance))+sum(PredictionError.*(PredictedObservedVariance\PredictionError),1)));
    % importance ratio 
    wtilde = lnw./wtilde ;
    % normalization 
    weights = wtilde/sum(wtilde);
    if (variance_update==1) && (neff(weights)<DynareOptions.particle.resampling.neff_threshold*sample_size)
        variance_update = 0 ;
    end
        % final resampling (advised)
    if second_resample==1 
        indx = index_resample(0,weights,DynareOptions);
        StateVectors = StateVectors(:,indx) ;
        if pruning
            StateVectors_ = StateVectors_(:,indx) ;
        end
        xparam = xparam(:,indx) ;
        weights = ones(1,number_of_particles)/number_of_particles ;
        mean_xparam(:,t) = mean(xparam,2);
        mat_var_cov = bsxfun(@minus,xparam,mean_xparam(:,t)) ;
        mat_var_cov = (mat_var_cov*mat_var_cov')/(number_of_particles-1) ;    
        std_xparam(:,t) = sqrt(diag(mat_var_cov)) ;
        for i=1:number_of_parameters
           temp = sortrows(xparam(i,:)') ;
           lb95_xparam(i,t) = temp(0.025*number_of_particles) ;
           median_xparam(i,t) = temp(0.5*number_of_particles) ;
           ub95_xparam(i,t) = temp(0.975*number_of_particles) ;
        end
    end
    if second_resample==0
        mean_xparam(:,t) = xparam*(weights') ;
        mat_var_cov = bsxfun(@minus,xparam,mean_xparam(:,t)) ;
        mat_var_cov = mat_var_cov*(bsxfun(@times,mat_var_cov,weights)') ;
        std_xparam(:,t) = sqrt(diag(mat_var_cov)) ;
        for i=1:number_of_parameters
           temp = sortrows([xparam(i,:)' weights'],1) ;
           cumulated_weights = cumsum(temp(:,2)) ;
           pass1=1 ;
           pass2=1 ;
           pass3=1 ;
           for j=1:number_of_particles
               if cumulated_weights(j)>0.025 && pass1==1 
                   lb95_xparam(i,t) = (temp(j-1,1)+temp(j,1))/2 ;
                   pass1 = 2 ;
               end
               if cumulated_weights(j)>0.5 && pass2==1
                   median_xparam(i,t) = (temp(j-1,1)+temp(j,1))/2 ;
                   pass2 = 2 ;
               end
               if cumulated_weights(j)>0.975 && pass3==1
                   ub95_xparam(i,t) = (temp(j-1,1)+temp(j,1))/2 ;
                   pass3 = 2 ;
               end
           end
        end
    end
    disp([lb95_xparam(:,t) mean_xparam(:,t) ub95_xparam(:,t)])
end
xparam = mean_xparam(:,sample_size) ; 
std_param = std_xparam(:,sample_size) ;
lb_95 = lb95_xparam(:,sample_size) ;
ub_95 = ub95_xparam(:,sample_size) ;
median_param = median_xparam(:,sample_size) ;

%% Plot parameters trajectory  
TeX = DynareOptions.TeX;

[nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_parameters);

if TeX
    fidTeX = fopen([Model.fname '_param_traj.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by online_auxiliary_filter.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

z = 1:1:sample_size ;

for plt = 1:nbplt,
    if TeX
        NAMES = [];
        TeXNAMES = [];
    end
    hh = dyn_figure(DynareOptions,'Name','Parameters Trajectories');
    for k=1:min(nstar,length(xparam)-(plt-1)*nstar)
        subplot(nr,nc,k)
        kk = (plt-1)*nstar+k;
        [name,texname] = get_the_name(kk,TeX,Model,EstimatedParameters,DynareOptions);
        if TeX
            if isempty(NAMES)
                NAMES = name;
                TeXNAMES = texname;
            else
                NAMES = char(NAMES,name);
                TeXNAMES = char(TeXNAMES,texname);
            end
        end
        y = [mean_xparam(kk,:)' median_xparam(kk,:)' lb95_xparam(kk,:)' ub95_xparam(kk,:)' xparam(kk)*ones(sample_size,1)] ;
        plot(z,y);
        hold on
        title(name,'interpreter','none')
        hold off
        axis tight
        drawnow
    end
    dyn_saveas(hh,[ Model.fname '_param_traj' int2str(plt) ],DynareOptions);
    if TeX
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:min(nstar,length(x)-(plt-1)*nstar)
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_ParamTraj%s}\n',Model.fname,int2str(plt));
        fprintf(fidTeX,'\\caption{Parameters trajectories.}');
        fprintf(fidTeX,'\\label{Fig:ParametersPlots:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
end

%% Plot Parameter Densities
number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourier Transform approximation.  
for plt = 1:nbplt,
    if TeX
        NAMES = [];
        TeXNAMES = [];
    end
    hh = dyn_figure(DynareOptions,'Name','Parameters Densities');
    for k=1:min(nstar,length(xparam)-(plt-1)*nstar)
        subplot(nr,nc,k)
        kk = (plt-1)*nstar+k;
        [name,texname] = get_the_name(kk,TeX,Model,EstimatedParameters,DynareOptions);
        if TeX
            if isempty(NAMES)
                NAMES = name;
                TeXNAMES = texname;
            else
                NAMES = char(NAMES,name);
                TeXNAMES = char(TeXNAMES,texname);
            end
        end
        optimal_bandwidth = mh_optimal_bandwidth(xparam(kk,:)',number_of_particles,bandwidth,kernel_function);
        [density(:,1),density(:,2)] = kernel_density_estimate(xparam(kk,:)',number_of_grid_points,...
                                                          number_of_particles,optimal_bandwidth,kernel_function);
        plot(density(:,1),density(:,2));
        hold on
        title(name,'interpreter','none')
        hold off
        axis tight
        drawnow
    end
    dyn_saveas(hh,[ Model.fname '_param_density' int2str(plt) ],DynareOptions);
    if TeX
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:min(nstar,length(x)-(plt-1)*nstar)
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_ParametersDensities%s}\n',Model.fname,int2str(plt));
        fprintf(fidTeX,'\\caption{ParametersDensities.}');
        fprintf(fidTeX,'\\label{Fig:ParametersDensities:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
end    
