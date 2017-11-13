function oo_ = realtime_shock_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% function oo_ = realtime_shock_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% Computes shocks contribution to a simulated trajectory. The fields set are
% oo_.realtime_shock_decomposition, oo_.conditional_shock_decomposition and oo_.realtime_forecast_shock_decomposition.
% Subfields are arrays n_var by nshock+2 by nperiods. The
% first nshock columns store the respective shock contributions, column n+1
% stores the role of the initial conditions, while column n+2 stores the
% value of the smoothed variables.  Both the variables and shocks are stored
% in the order of declaration, i.e. M_.endo_names and M_.exo_names, respectively.
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%    bayestopt_:  [structure]  describing the priors
%    estim_params_: [structure] characterizing parameters to be estimated
%
% OUTPUTS
%    oo_:         [structure]  Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2009-2017 Dynare Team
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

% indices of endogenous variables
if size(varlist,1) == 0
    varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
end

[i_var,nvar,index_uniques] = varlist_indices(varlist,M_.endo_names);
varlist=varlist(index_uniques,:);

% number of variables
endo_nbr = M_.endo_nbr;

% number of shocks
nshocks = M_.exo_nbr;

% parameter set
parameter_set = options_.parameter_set;
if isempty(parameter_set)
    if isfield(oo_,'posterior_mean')
        parameter_set = 'posterior_mean';
    elseif isfield(oo_,'mle_mode')
        parameter_set = 'mle_mode';
    elseif isfield(oo_,'posterior')
        parameter_set = 'posterior_mode';
    else
        error(['realtime_shock_decomposition: option parameter_set is not specified ' ...
               'and posterior mode is not available'])
    end
end

presample = max(1,options_.presample);
if isfield(options_.shock_decomp,'presample')
    presample = max(presample,options_.shock_decomp.presample);
end
% forecast_=0;
forecast_ = options_.shock_decomp.forecast;
forecast_params=0;
if forecast_ && isfield(options_.shock_decomp,'forecast_params')
    forecast_params = options_.shock_decomp.forecast_params;
end

% save_realtime=0;
save_realtime = options_.shock_decomp.save_realtime;
% array of time points in the range options_.presample+1:options_.nobs

zreal = zeros(endo_nbr,nshocks+2,options_.nobs+forecast_);
zcond = zeros(endo_nbr,nshocks+2,options_.nobs);
skipline()
skipline()
running_text = 'Realtime shock decomposition ';
newString=sprintf(running_text);
fprintf(['\b%s'],newString);

options_.selected_variables_only = 0; %make sure all variables are stored
options_.plot_priors=0;
init=1;
nobs = options_.nobs;

if forecast_ && any(forecast_params)
    M1=M_;
    M1.params = forecast_params;
    [junk1,junk2,junk3,junk4,junk5,junk6,oo1] = dynare_resolve(M1,options_,oo_);
    clear junk1 junk2 junk3 junk4 junk5 junk6
end

for j=presample+1:nobs
    %    evalin('base',['options_.nobs=' int2str(j) ';'])
    options_.nobs=j;
    [oo, M_, junk2, junk3, Smoothed_Variables_deviation_from_mean] = evaluate_smoother(parameter_set,varlist,M_,oo_,options_,bayestopt_,estim_params_);

    % reduced form
    dr = oo.dr;

    % data reordering
    order_var = dr.order_var;
    inv_order_var = dr.inv_order_var;


    % coefficients
    A = dr.ghx;
    B = dr.ghu;

    if forecast_
        if any(forecast_params)
            Af = oo1.dr.ghx;
            Bf = oo1.dr.ghu;
        else
            Af = A;
            Bf = B;
        end
    end

    % initialization
    gend = size(oo.SmoothedShocks.(deblank(M_.exo_names(1,:))),1);
    epsilon=NaN(nshocks,gend);
    for i=1:nshocks
        epsilon(i,:) = oo.SmoothedShocks.(deblank(M_.exo_names(i,:)));
    end
    epsilon=[epsilon zeros(nshocks,forecast_)];

    z = zeros(endo_nbr,nshocks+2,gend+forecast_);

    z(:,end,1:gend) = Smoothed_Variables_deviation_from_mean;

    maximum_lag = M_.maximum_lag;

    k2 = dr.kstate(find(dr.kstate(:,2) <= maximum_lag+1),[1 2]);
    i_state = order_var(k2(:,1))+(min(i,maximum_lag)+1-k2(:,2))*M_.endo_nbr;
    for i=1:gend+forecast_
        if i > 1 && i <= maximum_lag+1
            lags = min(i-1,maximum_lag):-1:1;
        end

        if i > 1
            tempx = permute(z(:,1:nshocks,lags),[1 3 2]);
            m = min(i-1,maximum_lag);
            tempx = [reshape(tempx,endo_nbr*m,nshocks); zeros(endo_nbr*(maximum_lag-i+1),nshocks)];
            if i > gend
                z(:,nshocks+2,i) = Af(inv_order_var,:)*z(i_state,nshocks+2,lags);
                %             z(:,nshocks+2,i) = A(inv_order_var,:)*permute(z(i_state,nshocks+2,lags),[1 3 2]);
                z(:,1:nshocks,i) = Af(inv_order_var,:)*tempx(i_state,:);
            else
                z(:,1:nshocks,i) = A(inv_order_var,:)*tempx(i_state,:);
            end
            lags = lags+1;
            z(:,1:nshocks,i) = z(:,1:nshocks,i) + B(inv_order_var,:).*repmat(epsilon(:,i)',endo_nbr,1);
        end

        %         z(:,1:nshocks,i) = z(:,1:nshocks,i) + B(inv_order_var,:).*repmat(epsilon(:,i)',endo_nbr,1);
        z(:,nshocks+1,i) = z(:,nshocks+2,i) - sum(z(:,1:nshocks,i),2);
    end

    %% conditional shock decomp 1 step ahead
    z1 = zeros(endo_nbr,nshocks+2);
    z1(:,end) = Smoothed_Variables_deviation_from_mean(:,gend);
    for i=gend

        z1(:,1:nshocks) = z1(:,1:nshocks) + B(inv_order_var,:).*repmat(epsilon(:,i)',endo_nbr,1);
        z1(:,nshocks+1) = z1(:,nshocks+2) - sum(z1(:,1:nshocks),2);
    end
    %%

    %% conditional shock decomp k step ahead
    if forecast_ && forecast_<j
        zn = zeros(endo_nbr,nshocks+2,forecast_+1);
        zn(:,end,1:forecast_+1) = Smoothed_Variables_deviation_from_mean(:,gend-forecast_:gend);
        for i=1:forecast_+1
            if i > 1 && i <= maximum_lag+1
                lags = min(i-1,maximum_lag):-1:1;
            end

            if i > 1
                tempx = permute(zn(:,1:nshocks,lags),[1 3 2]);
                m = min(i-1,maximum_lag);
                tempx = [reshape(tempx,endo_nbr*m,nshocks); zeros(endo_nbr*(maximum_lag-i+1-1),nshocks)];
                zn(:,1:nshocks,i) = A(inv_order_var,:)*tempx(i_state,:);
                lags = lags+1;
                zn(:,1:nshocks,i) = zn(:,1:nshocks,i) + B(inv_order_var,:).*repmat(epsilon(:,i+gend-forecast_-1)',endo_nbr,1);
            end

            %             zn(:,1:nshocks,i) = zn(:,1:nshocks,i) + B(inv_order_var,:).*repmat(epsilon(:,i+gend-forecast_-1)',endo_nbr,1);
            zn(:,nshocks+1,i) = zn(:,nshocks+2,i) - sum(zn(:,1:nshocks,i),2);
        end
        oo_.conditional_shock_decomposition.(['time_' int2str(j-forecast_)])=zn;
    end
    %%

    if init
        zreal(:,:,1:j) = z(:,:,1:j);
    else
        zreal(:,:,j) = z(:,:,gend);
    end
    zcond(:,:,j) = z1;
    if ismember(j,save_realtime)
        oo_.realtime_shock_decomposition.(['time_' int2str(j)])=z;
    end

    if forecast_
        zfrcst(:,:,j+1) = z(:,:,gend+1);
        oo_.realtime_forecast_shock_decomposition.(['time_' int2str(j)])=z(:,:,gend:end);
        if j>forecast_+presample
            %% realtime conditional shock decomp k step ahead
            oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-forecast_)]) = ...
                zreal(:,:,j-forecast_:j) - ...
                oo_.realtime_forecast_shock_decomposition.(['time_' int2str(j-forecast_)]);
            oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-forecast_)])(:,end-1,:) = ...
                oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-forecast_)])(:,end-1,:) + ...
                oo_.realtime_forecast_shock_decomposition.(['time_' int2str(j-forecast_)])(:,end,:);
            oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-forecast_)])(:,end,:) = ...
                zreal(:,end,j-forecast_:j);

            if j==nobs
                for my_forecast_=(forecast_-1):-1:1
                    oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-my_forecast_)]) = ...
                        zreal(:,:,j-my_forecast_:j) - ...
                        oo_.realtime_forecast_shock_decomposition.(['time_' int2str(j-my_forecast_)])(:,:,1:my_forecast_+1);
                    oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-my_forecast_)])(:,end-1,:) = ...
                        oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-my_forecast_)])(:,end-1,:) + ...
                        oo_.realtime_forecast_shock_decomposition.(['time_' int2str(j-my_forecast_)])(:,end,1:my_forecast_+1);
                    oo_.realtime_conditional_shock_decomposition.(['time_' int2str(j-my_forecast_)])(:,end,:) = ...
                        zreal(:,end,j-my_forecast_:j);
                end
            end

        end
    end

    prctdone=(j-presample)/(nobs-presample);
    if isoctave
        printf([running_text,' %3.f%% done\r'], prctdone*100);
    else
        s0=repmat('\b',1,length(newString));
        newString=sprintf([running_text,' %3.1f%% done'], prctdone*100);
        fprintf([s0,'%s'],newString);
    end
    init=0;
end
oo_.realtime_shock_decomposition.pool = zreal;
oo_.conditional_shock_decomposition.pool = zcond;
if forecast_
    oo_.realtime_forecast_shock_decomposition.pool = zfrcst;
end

skipline()
