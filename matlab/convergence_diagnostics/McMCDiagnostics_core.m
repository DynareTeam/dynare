function myoutput = McMCDiagnostics_core(myinputs,fpar,npar,whoiam, ThisMatlab)
% function myoutput = McMCDiagnostics_core(myinputs,fpar,npar,whoiam, ThisMatlab)
% Computes the Brooks/Gelman (1998) convergence diagnostics, both the
% parameteric and the non-parameteric versions
% 
% PARALLEL CONTEXT
% Core functionality for MCMC Diagnostics, which can be parallelized.
% See also the comment in posterior_sampler_core.m funtion.
% 
% 
% INPUTS
%   See See the comment in posterior_sampler_core.m funtion.

% OUTPUTS
% o myoutput  [struc]
%   Contains:
%       - UDIAG [by 6] double   1st column: length of total sequence interval
%                               2nd column: sum of length of within sequence intervals; used to compute mean length of within sequence intervals
%                               3nd column: within sequence variance
%                               4nd column: sum of within sequence variances; used to compute mean within sequence variances
%                               5nd column: within sequence kurtosis
%                               6nd column: sum of within sequence kurtoses; used to compute mean within sequence kurtoses
%               Averaging to compute mean moments is done in
%               McMCDiagnostics
%
% ALGORITHM
%   Computes part of the convergence diagnostics, the rest is computed in McMCDiagnostics.m .
%   The methodology and terminology is based on: Brooks/Gelman (1998): General 
%   Methods for Monitoring Convergence of Iterative Simulations, Journal of Computational 
%   and Graphical Statistics, Volume 7, Number 4, Pages 434-455
% 
%
% SPECIAL REQUIREMENTS.
%   None.

% Copyright (C) 2006-2016 Dynare Team
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

if nargin<4,
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

MetropolisFolder=myinputs.MetropolisFolder;%myinputs.MetropolisFolder;
nblck=myinputs.nblck;
NumberOfMcFilesPerBlock=myinputs.NumberOfMcFilesPerBlock;
Origin=myinputs.Origin;
StepSize=myinputs.StepSize;
mh_drop=myinputs.mh_drop;
NumberOfDraws=myinputs.NumberOfDraws;
NumberOfLines=myinputs.NumberOfLines;
time=myinputs.time;
M_=myinputs.M_;

if whoiam
    Parallel=myinputs.Parallel;
end
if ~exist('MetropolisFolder'),
    MetropolisFolder = CheckPath('metropolis',M_.dname);
end

ALPHA = 0.2; % percentile for non-parametric statistic
tmp = zeros(NumberOfDraws*nblck,3);
UDIAG = zeros(NumberOfLines,6,npar-fpar+1);

if whoiam
    waitbarString = ['Please wait... McMCDiagnostics (' int2str(fpar) 'of' int2str(npar) ')...'];
    if Parallel(ThisMatlab).Local,
        waitbarTitle=['Local '];
    else
        waitbarTitle=[Parallel(ThisMatlab).ComputerName];
    end
    fMessageStatus(0,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
end
for j=fpar:npar,
    if isoctave
        if (whoiam==0),
            printf('    Parameter %d...  ',j);
        end
    else
        fprintf('    Parameter %d...  ',j);
    end
    for b = 1:nblck %load draws from different chains into 1 matrix
        startline = 0;
        for n = 1:NumberOfMcFilesPerBlock
            load([MetropolisFolder '/' M_.fname '_mh',int2str(n),'_blck' int2str(b) '.mat'],'x2');
            nx2 = size(x2,1);
            tmp((b-1)*NumberOfDraws+startline+(1:nx2),1) = x2(:,j);
            startline = startline + nx2;
        end
    end
    tmp(:,2) = kron(transpose(1:nblck),ones(NumberOfDraws,1)); %add info about chain associated with draw into 2nd column
    tmp(:,3) = kron(ones(nblck,1),time'); %add timeline for draws to third column
    tmp = sortrows(tmp,1); %sort draws according to size for non-parametric percentile computation
    window_iter   = 0;
    for iter  = Origin:StepSize:NumberOfDraws %begin of window
        window_iter = window_iter+1;
        linea = ceil(mh_drop*iter);         %compute first non-discarded draw; drops fraction of sample at each iteration for computational efficiency, see Brooks/Gelman (1998), p.438
        n     = iter-linea+1;               %number of draws from each block in current batch
        cinf  = round(n*ALPHA/2);           %lower bound for alpha percentile of within series
        csup  = round(n*(1-ALPHA/2));       %upper bound for alpha percentile of within series
        CINF  = round(nblck*n*ALPHA/2);     %lower bound for alpha percentile of pooled series
        CSUP  = round(nblck*n*(1-ALPHA/2)); %upper bound for alpha percentile of pooled series
        temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);    %extract pooled draws in current batch
        UDIAG(window_iter,1,j-fpar+1) = temp(CSUP,1)-temp(CINF,1);      %length of total sequence interval
        pooled_mean = mean(temp(:,1));      % Pooled mean.
        UDIAG(window_iter,3,j-fpar+1) = sum((temp(:,1)-pooled_mean).^2)/(nblck*n-1);        %within sequence variance
        UDIAG(window_iter,5,j-fpar+1) = sum(abs(temp(:,1)-pooled_mean).^3)/(nblck*n-1);     %within sequence third moment
        for i=1:nblck
            pmet = temp(find(temp(:,2)==i));
            UDIAG(window_iter,2,j-fpar+1) = UDIAG(window_iter,2,j-fpar+1) + pmet(csup,1)-pmet(cinf,1); %sum of length of within sequence intervals; used to compute mean length of within sequence intervals
            within_mean = mean(pmet,1); %% Within mean in current chain.
            UDIAG(window_iter,4,j-fpar+1) = UDIAG(window_iter,4,j-fpar+1) + sum((pmet(:,1)-within_mean).^2)/(n-1); %sum of within sequence variances; used to compute mean within sequence variances
            UDIAG(window_iter,6,j-fpar+1) = UDIAG(window_iter,6,j-fpar+1) + sum(abs(pmet(:,1)-within_mean).^3)/(n-1); %sum of within sequence kurtoses; used to compute mean within sequence kurtoses
        end
    end
    if isoctave
        if (whoiam==0),
            printf('Done! \n');
        end
    else
        fprintf('Done! \n');
    end
    if whoiam,
        waitbarString = [ 'Parameter ' int2str(j) '/' int2str(npar) ' done.'];
        fMessageStatus((j-fpar+1)/(npar-fpar+1),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab))
    end
end

myoutput.UDIAG = UDIAG;