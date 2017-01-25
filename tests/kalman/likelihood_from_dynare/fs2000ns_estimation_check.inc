%%get mode
estimation(diffuse_filter,kalman_algo=3,mode_compute=4,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_0=oo_.likelihood_at_initial_parameters;

%%Diffuse Multivariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=3,diffuse_filter,mode_file=@{mode_file_name},mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_2=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,1)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,1)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,1)=cell2mat(struct2cell(oo_.SmoothedVariables));

%%Diffuse univariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=4,diffuse_filter,mode_file=@{mode_file_name},mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_4=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,2)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,2)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,2)=cell2mat(struct2cell(oo_.SmoothedVariables));

if max(max(abs(SmoothedMeasurementErrors-repmat(SmoothedMeasurementErrors(:,:,1),[1,1,2]))))>1e-8
    error('SmoothedMeasurementErrors do not match')
end

if max(max(abs(SmoothedShocks-repmat(SmoothedShocks(:,:,1),[1,1,2]))))>1e-8
    error('SmoothedShocks do not match')
end

if isoctave
    tol_crit=5e-6;
else
    tol_crit=1e-8;    
end

if max(max(abs(SmoothedVariables-repmat(SmoothedVariables(:,:,1),[1,1,2]))))>tol_crit
    error('SmoothedVariables do not match')
end

if max(abs(fval_algo_2-fval_algo_4))>1e-6
    fprintf('\nLikelihood difference: %4.3f\n',max(abs(fval_algo_2-fval_algo_4)))
    error('Likelihoods do not match')
end
