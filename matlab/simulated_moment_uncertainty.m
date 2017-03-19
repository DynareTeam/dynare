function [cmm, mm] = simulated_moment_uncertainty(indx, periods, replic,options_,M_,oo_)
% function [cmm, mm] = simulated_moment_uncertainty(indx, periods, replic,options_,M_,oo_)
% Compute the uncertainty around simulated moments
% Inputs
%   - indx      [n_moments by 1]  index vector of moments
%   - periods   [scalar]    number of simulation periods
%   - replic    [scalar]    number of simulation replications
%   - options_  Dynare options structure
%   - M_        Dynare Model structure
%   - oo_       Dynare results structure
% Outputs:
%   - cmm:      [n_moments by n_moments] covariance matrix of simulated moments
%   - mm:       [n_moments by replic] matrix of moments
%  
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

mm=zeros(length(indx),replic);

disp('Evaluating simulated moment uncertainty ... please wait')
disp(['Doing ',int2str(replic),' replicas of length ',int2str(periods),' periods.'])

h = dyn_waitbar(0,'Simulated moment uncertainty ...');
%Do check whether simulation is possible
if options_.periods == 0
    error('simulated_moment_uncertainty: Periods must be bigger than 0')
end
if options_.periods <= options_.drop
    error('simulated_moment_uncertainty: The horizon of simulation is shorter than the number of observations to be dropped. Either increase options_.periods or decrease options_.drop.')
end

%locally set options
options_.TeX=0;
options_.noprint = 1;
options_.order = 1;
options_.periods = periods;

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end


if M_.exo_nbr > 0
    oo_.exo_simul= ones(max(options_.periods,1) + M_.maximum_lag + M_.maximum_lead,1) * oo_.exo_steady_state';
end

oo_.dr=set_state_space(oo_.dr,M_,options_);


if options_.logged_steady_state %if steady state was previously logged, undo this
    oo_.dr.ys=exp(oo_.dr.ys);
    oo_.steady_state=exp(oo_.steady_state);
    options_.logged_steady_state=0;
    logged_steady_state_indicator=1;
    evalin('base','options_.logged_steady_state=0;')
else 
    logged_steady_state_indicator=0;
end

[dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
oo_.dr=dr;
if info(1)
    fprintf('\nsimulated_moment_uncertainty: model could not be solved')
    print_info(info,0,options_);
end

%set starting point of simulations
if isempty(M_.endo_histval)
    if options_.loglinear
        y0 = log(oo_.dr.ys);
    else
        y0 = oo_.dr.ys;
    end
else
    if options_.loglinear
        y0 = log_variable(1:M_.endo_nbr,M_.endo_histval,M_);
    else
        y0 = M_.endo_histval;
    end
end


for j=1:replic;
    [ys, oo_] = simult(y0,oo_.dr,M_,options_,oo_);%do simulation
    oo_=disp_moments(ys,char(options_.varobs),M_,options_,oo_); %get moments
    dum=[oo_.mean; dyn_vech(oo_.var)];
    sd = sqrt(diag(oo_.var));
    for i=1:options_.ar;
        dum=[dum; vec(oo_.autocorr{i}.*(sd*sd'))];
    end
    mm(:,j)=dum(indx);
    dyn_waitbar(j/replic,h,['Simulated moment uncertainty. Replic  ',int2str(j),'/',int2str(replic)])
end;
dyn_waitbar_close(h);

if logged_steady_state_indicator
    evalin('base','options_.logged_steady_state=1;') %reset base workspace option to conform to base oo_
end
cmm = cov(mm');
disp('Simulated moment uncertainty ... done!')
