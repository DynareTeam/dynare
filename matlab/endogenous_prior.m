function [lnpriormom] = endogenous_prior(data,Pstar,BayesInfo,H)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
%    data           [double]     n*T vector of data observations 
%    Pstar          [double]     k*k matrix of
%    BayesInfo      [structure]   
%
% OUTPUTS
%    lnpriormom     [double]     scalar of log endogenous prior value

% Code to implement notes on endogenous priors by Lawrence Christiano,
% specified in the appendix of:
% ’Introducing Financial Frictions and Unemployment into a Small Open Economy Model’
% by Lawrence J. Christiano, Mathias Trabandt and Karl Walentin (2011), Journal of Economic Dynamics and Control      
% this is the 'mother' of the priors on the model parameters.
% the priors include a metric across some choosen moments of the (supposedly
% pre-sample) data.
% *** Implemented file for variances, but in principle any moment
% *** could be matched
% As a default, the prior second moments are computed from the same sample
% used to find the posterior mode. This could be changed by making the
% appropriate adjustment to the following code.


% Copyright (C) 2011 Lawrence J. Christiano, Mathias Trabandt and Karl Walentin 
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

Y=data';
[Tsamp,n]=size(Y);    % sample length and number of matched moments (here set equal to nr of observables)

hmat=zeros(n,Tsamp);
Ydemean=zeros(Tsamp,n);
C0=zeros(n,n);
C1=zeros(n,n);
C2=zeros(n,n);

for j=1:n
    Ydemean(:,j)=Y(:,j)-mean(Y(:,j));
end
Fhat=diag(Ydemean'*Ydemean)/Tsamp;

% we need ht, where t=1,...,T
for t=1:Tsamp
    hmat(:,t)=diag(Ydemean(t,:)'*Ydemean(t,:))-Fhat;
end

% To calculate Shat we need C0, C1 and C2
for t=1:Tsamp
    C0=C0+1/Tsamp*hmat(:,t)*hmat(:,t)';
end

for t=2:Tsamp
    C1=C1+1/(Tsamp-1)*hmat(:,t)*hmat(:,t-1)';
end

for t=3:Tsamp
    C2=C2+1/(Tsamp-2)*hmat(:,t)*hmat(:,t-2)';
end

% Finally, we have the sampling uncertainty measure Shat:
Shat=C0 +(1-1/(2+1))*(C1+C1')...
    +(1-2/(2+1))*(C2+C2');

% Model variances below:
mf=BayesInfo.mf1;
II=eye(size(Pstar,2));
Z=II(mf,:);
% This is Ftheta, variance of model variables, given param vector theta:
Ftheta=diag(Z*Pstar(:,mf)+H);
% below commented out line is for Del Negro Schorfheide style priors:
%     lnpriormom=-.5*n*TT*log(2*pi)-.5*TT*log(det(sigma))-.5*TT*trace(inv(sigma)*(gamyy-2*phi'*gamxy+phi'*gamxx*phi));
lnpriormom=.5*n*log(Tsamp/(2*pi))-.5*log(det(Shat))-.5*Tsamp*(Fhat-Ftheta)'/Shat*(Fhat-Ftheta);

