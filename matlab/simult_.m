function y_=simult_(y0,dr,ex_,iorder)
% Simulates the model using a perturbation approach, given the path for the exogenous variables and the
% decision rules.
%
% INPUTS
%    y0       [double]   n*1 vector, initial value (n is the number of declared endogenous variables plus the number 
%                        of auxilliary variables for lags and leads)
%    dr       [struct]   matlab's structure where the reduced form solution of the model is stored.
%    ex_      [double]   T*q matrix of innovations.
%    iorder   [integer]  order of the taylor approximation.
%
% OUTPUTS
%    y_       [double]   n*(T+1) time series for the endogenous variables.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2012 Dynare Team
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

global M_ options_

iter = size(ex_,1);
endo_nbr = M_.endo_nbr;
exo_nbr = M_.exo_nbr;

y_ = zeros(size(y0,1),iter+M_.maximum_lag);
y_(:,1) = y0;

% stoch_simul sets k_order_solver=1 if order=3, but does so only locally, so we
% have to do it here also
if options_.order == 3
    options_.k_order_solver = 1;
end

if ~options_.k_order_solver || (options_.k_order_solver && options_.pruning) %if k_order_pert is not used or if we do not use Dynare++ with k_order_pert
    if iorder==1
        y_(:,1) = y_(:,1)-dr.ys;
    end
end

if options_.k_order_solver && ~options_.pruning % Call dynare++ routines.
    ex_ = [zeros(1,exo_nbr); ex_];
    switch options_.order
      case 1
        [err, y_] = dynare_simul_(1,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                                  y_(dr.order_var,1),ex_',M_.Sigma_e,options_.DynareRandomStreams.seed,dr.ys(dr.order_var),...
                                  zeros(endo_nbr,1),dr.g_1);
      case 2
        [err, y_] = dynare_simul_(2,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                                  y_(dr.order_var,1),ex_',M_.Sigma_e,options_.DynareRandomStreams.seed,dr.ys(dr.order_var),dr.g_0, ...
                                  dr.g_1,dr.g_2);
      case 3
        [err, y_] = dynare_simul_(3,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                                  y_(dr.order_var,1),ex_',M_.Sigma_e,options_.DynareRandomStreams.seed,dr.ys(dr.order_var),dr.g_0, ...
                                  dr.g_1,dr.g_2,dr.g_3);
      otherwise
        error(['order = ' int2str(order) ' isn''t supported'])
    end
    mexErrCheck('dynare_simul_', err);
    y_(dr.order_var,:) = y_;
else
    if options_.block
        if M_.maximum_lag > 0
            k2 = dr.state_var;
        else
            k2 = [];
        end;
        order_var = 1:endo_nbr;
        dr.order_var = order_var;
    else
        k2 = dr.kstate(find(dr.kstate(:,2) <= M_.maximum_lag+1),[1 2]);
        k2 = k2(:,1)+(M_.maximum_lag+1-k2(:,2))*endo_nbr;
        order_var = dr.order_var;
    end;
    
    switch iorder
      case 1
        if isempty(dr.ghu)% For (linearized) deterministic models.
            for i = 2:iter+M_.maximum_lag
                yhat = y_(order_var(k2),i-1);
                y_(order_var,i) = dr.ghx*yhat;
            end
        elseif isempty(dr.ghx)% For (linearized) purely forward variables (no state variables).
            y_(dr.order_var,:) = dr.ghu*transpose(ex_);
        else
            epsilon = dr.ghu*transpose(ex_);
            for i = 2:iter+M_.maximum_lag
                yhat = y_(order_var(k2),i-1);
                y_(order_var,i) = dr.ghx*yhat + epsilon(:,i-1);
            end
        end
        y_ = bsxfun(@plus,y_,dr.ys);
      case 2
        constant = dr.ys(order_var)+.5*dr.ghs2;
        if options_.pruning
            y__ = y0;
            for i = 2:iter+M_.maximum_lag
                yhat1 = y__(order_var(k2))-dr.ys(order_var(k2));
                yhat2 = y_(order_var(k2),i-1)-dr.ys(order_var(k2));
                epsilon = ex_(i-1,:)';
                [abcOut1, err] = A_times_B_kronecker_C(.5*dr.ghxx,yhat1,options_.threads.kronecker.A_times_B_kronecker_C);
                mexErrCheck('A_times_B_kronecker_C', err);
                [abcOut2, err] = A_times_B_kronecker_C(.5*dr.ghuu,epsilon,options_.threads.kronecker.A_times_B_kronecker_C);
                mexErrCheck('A_times_B_kronecker_C', err);
                [abcOut3, err] = A_times_B_kronecker_C(dr.ghxu,yhat1,epsilon,options_.threads.kronecker.A_times_B_kronecker_C);
                mexErrCheck('A_times_B_kronecker_C', err);
                y_(order_var,i) = constant + dr.ghx*yhat2 + dr.ghu*epsilon ...
                    + abcOut1 + abcOut2 + abcOut3;
                y__(order_var) = dr.ys(order_var) + dr.ghx*yhat1 + dr.ghu*epsilon;
            end
        else
            for i = 2:iter+M_.maximum_lag
                yhat = y_(order_var(k2),i-1)-dr.ys(order_var(k2));
                epsilon = ex_(i-1,:)';
                [abcOut1, err] = A_times_B_kronecker_C(.5*dr.ghxx,yhat,options_.threads.kronecker.A_times_B_kronecker_C);
                mexErrCheck('A_times_B_kronecker_C', err);
                [abcOut2, err] = A_times_B_kronecker_C(.5*dr.ghuu,epsilon,options_.threads.kronecker.A_times_B_kronecker_C);
                mexErrCheck('A_times_B_kronecker_C', err);
                [abcOut3, err] = A_times_B_kronecker_C(dr.ghxu,yhat,epsilon,options_.threads.kronecker.A_times_B_kronecker_C);
                mexErrCheck('A_times_B_kronecker_C', err);
                y_(dr.order_var,i) = constant + dr.ghx*yhat + dr.ghu*epsilon ...
                    + abcOut1 + abcOut2 + abcOut3;
            end
         end
      case 3
        % only with pruning
        % the third moments of the shocks are assumed null. We don't have
        % an interface for specifying them
        ghx = dr.ghx;
        ghu = dr.ghu;
        ghxx = dr.ghxx;
        ghxu = dr.ghxu;
        ghuu = dr.ghuu;
        ghs2 = dr.ghs2;
        ghxxx = dr.ghxxx;
        ghxxu = dr.ghxxu;
        ghxuu = dr.ghxuu;
        ghuuu = dr.ghuuu;
        ghxss = dr.ghxss;
        ghuss = dr.ghuss;
        threads = options_.threads.kronecker.A_times_B_kronecker_C;
        nspred = M_.nspred;
        ipred = M_.nstatic+(1:nspred);
        %construction follows Andreasen et al (2013), Technical
        %Appendix, Formulas (65) and (66)
        %split into first, second, and third order terms
        yhat1 = y0(order_var(k2))-dr.ys(order_var(k2));
        yhat2 = zeros(nspred,1);
        yhat3 = zeros(nspred,1);
        for i=2:iter+M_.maximum_lag
            u = ex_(i-1,:)';
            %construct terms of order 2 from second order part, based
            %on linear component yhat1
            [gyy, err] = A_times_B_kronecker_C(ghxx,yhat1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [guu, err] = A_times_B_kronecker_C(ghuu,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gyu, err] = A_times_B_kronecker_C(ghxu,yhat1,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            %construct terms of order 3 from second order part, based
            %on order 2 component yhat2
            [gyy12, err] = A_times_B_kronecker_C(ghxx,yhat1,yhat2,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gy2u, err] = A_times_B_kronecker_C(ghxu,yhat2,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            %construct terms of order 3, all based on first order component yhat1              
            y2a = kron(yhat1,yhat1);
            [gyyy, err] = A_times_B_kronecker_C(ghxxx,y2a,yhat1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            u2a = kron(u,u);
            [guuu, err] = A_times_B_kronecker_C(ghuuu,u2a,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            yu = kron(yhat1,u);
            [gyyu, err] = A_times_B_kronecker_C(ghxxu,yhat1,yu,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gyuu, err] = A_times_B_kronecker_C(ghxuu,yu,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            %add all terms of order 3, linear component based on third
            %order yhat3
            yhat3 = ghx*yhat3 +gyy12 ... % prefactor is 1/2*2=1, see (65) Appendix Andreasen et al.
                    + gy2u ... % prefactor is 1/2*2=1, see (65) Appendix Andreasen et al.
                    + 1/6*(gyyy + guuu + 3*(gyyu + gyuu +  ghxss*yhat1 + ghuss*u)); %note: s is treated as variable, thus xss and uss are third order
            yhat2 = ghx*yhat2 + 1/2*(gyy + guu + 2*gyu + ghs2);
            yhat1 = ghx*yhat1 + ghu*u;
            y_(order_var,i) = dr.ys(order_var)+yhat1 + yhat2 + yhat3; %combine terms again
            yhat1 = yhat1(ipred);
            yhat2 = yhat2(ipred);
            yhat3 = yhat3(ipred);
        end   
    end
end
