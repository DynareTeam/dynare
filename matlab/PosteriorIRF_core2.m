function myoutput=PosteriorIRF_core2(myinputs,fpar,npar,whoiam, ThisMatlab)
% Generates the Posterior IRFs plot from the IRFs generated in
% PosteriorIRF_core1
% PARALLEL CONTEXT
% Perform in parallel execution a portion of the PosteriorIRF.m code.
% See also the comment in random_walk_metropolis_hastings_core.m funtion.
%
% INPUTS 
%   See the comment in random_walk_metropolis_hastings_core.m funtion.
%
% OUTPUTS
% o myoutput  [struc]
%  Contained:
%  OutputFileName (i.e. the figures without the file .txt).
%
% ALGORITHM 
%   Portion of PosteriorIRF.m function code. Specifically the last 'for' cycle.       
%
% SPECIAL REQUIREMENTS.
%   None.
%
% Copyright (C) 2006-2013 Dynare Team
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

global options_  M_ 

if nargin<4,
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

Check=options_.TeX;
if (Check)
    varlist_TeX=myinputs.varlist_TeX;
end

nvar=myinputs.nvar;
MeanIRF=myinputs.MeanIRF;
tit=myinputs.tit;
nn=myinputs.nn;
MAX_nirfs_dsgevar=myinputs.MAX_nirfs_dsgevar;
HPDIRF=myinputs.HPDIRF;
if options_.dsge_var
    HPDIRFdsgevar=myinputs.HPDIRFdsgevar;
    MeanIRFdsgevar=myinputs.MeanIRFdsgevar;
end

varlist=myinputs.varlist;
MaxNumberOfPlotPerFigure=myinputs.MaxNumberOfPlotPerFigure;

% Necessary only for remote computing!
if whoiam
    Parallel=myinputs.Parallel;
end

% To save the figures where the function is computed!

DirectoryName = CheckPath('Output',M_.dname);

RemoteFlag = 0;
if whoiam,
    if Parallel(ThisMatlab).Local==0,
        RemoteFlag =1;
    end
    prct0={0,whoiam,Parallel(ThisMatlab)};
    dyn_waitbar(prct0,'PosteriorIRF Plots ...');
end

OutputFileName={};

subplotnum = 0;
for i=fpar:npar,
    figunumber = 0;
    
    for j=1:nvar
        if max(abs(MeanIRF(:,j,i))) > options_.impulse_responses.plot_threshold
            subplotnum = subplotnum+1;
            if subplotnum == 1 && options_.relative_irf
                hh = dyn_figure(options_,'Name',['Relative response to orthogonalized shock to ' tit(i,:)]);
            elseif subplotnum == 1 && ~options_.relative_irf
                hh = dyn_figure(options_,'Name',['Orthogonalized shock to ' tit(i,:)]);
            end
            
            set(0,'CurrentFigure',hh)
            subplot(nn,nn,subplotnum);
            if ~MAX_nirfs_dsgevar
                h1 = area(1:options_.irf,HPDIRF(:,2,j,i));
                set(h1,'FaceColor',[.9 .9 .9]);
                set(h1,'BaseValue',min(HPDIRF(:,1,j,i)));
                hold on
                h2 = area(1:options_.irf,HPDIRF(:,1,j,i),'FaceColor',[1 1 1],'BaseValue',min(HPDIRF(:,1,j,i)));
                set(h2,'FaceColor',[1 1 1]);
                set(h2,'BaseValue',min(HPDIRF(:,1,j,i)));
                plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',3)
                % plot([1 options_.irf],[0 0],'-r','linewidth',0.5);          
                box on
                axis tight
                xlim([1 options_.irf]);
                hold off
            else    
                h1 = area(1:options_.irf,HPDIRF(:,2,j,i));
                set(h1,'FaceColor',[.9 .9 .9]);
                set(h1,'BaseValue',min([min(HPDIRF(:,1,j,i)),min(HPDIRFdsgevar(:,1,j,i))]));
                hold on;
                h2 = area(1:options_.irf,HPDIRF(:,1,j,i));
                set(h2,'FaceColor',[1 1 1]);
                set(h2,'BaseValue',min([min(HPDIRF(:,1,j,i)),min(HPDIRFdsgevar(:,1,j,i))]));
                plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',3)
                % plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
                plot(1:options_.irf,MeanIRFdsgevar(:,j,i),'--k','linewidth',2)
                plot(1:options_.irf,HPDIRFdsgevar(:,1,j,i),'--k','linewidth',1)
                plot(1:options_.irf,HPDIRFdsgevar(:,2,j,i),'--k','linewidth',1)
                box on
                axis tight
                xlim([1 options_.irf]);
                hold off
            end
            name = deblank(varlist(j,:));
            title(name,'Interpreter','none')
        else
            if options_.debug
                fprintf('POSTERIOR_IRF: The IRF of %s to %s is smaller than the irf_plot_threshold of %4.3f and will not be displayed.\n',deblank(varlist(j,:)),tit(i,:),options_.impulse_responses.plot_threshold)
            end                
        end
        
        if subplotnum == MaxNumberOfPlotPerFigure || (j == nvar  && subplotnum> 0)
            figunumber = figunumber+1;
            dyn_saveas(hh,[DirectoryName '/'  M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber)],options_);
            if RemoteFlag==1,
                OutputFileName = [OutputFileName; {[DirectoryName,filesep], [M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber) '.*']}];
            end
            subplotnum = 0;
        end
    end% loop over selected endo_var
    if whoiam,
        fprintf('Done! \n');
        waitbarString = [ 'Exog. shocks ' int2str(i) '/' int2str(npar) ' done.'];
%         fMessageStatus((i-fpar+1)/(npar-fpar+1),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab));
        dyn_waitbar((i-fpar+1)/(npar-fpar+1),[],waitbarString);
    end
end% loop over exo_var  



myoutput.OutputFileName = OutputFileName;

