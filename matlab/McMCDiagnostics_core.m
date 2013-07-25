function myoutput = McMCDiagnostics_core(myinputs,fpar,npar,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% Core functionality for MCMC Diagnostics, which can be parallelized.
% See also the comment in random_walk_metropolis_hastings_core.m funtion.


% INPUTS
%   See See the comment in random_walk_metropolis_hastings_core.m funtion.

% OUTPUTS
% o myoutput  [struc]
%  Contained UDIAG.
%
% ALGORITHM
%   Portion of McMCDiagnostics.m function.
%
% SPECIAL REQUIREMENTS.
%   None.

% Copyright (C) 2006-2011 Dynare Team
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

%%
% % Unpack Global Variables, they are already in global space when running
% % single computation
% % 
try                                                 % would be catched in single computation ('no such fieldname')
    globalVars = fieldnames(myinputs.global);       % packed by masterParallel2      
    for j=1:length(globalVars),
        eval(['global ',globalVars{j},';'])
        fieldname=globalVars{j};
        value=myinputs.global.(fieldname);
        eval([fieldname '=value;'])
        
        evalin('base',['global ', globalVars{j},';']) % put also into base workspace
        assignin('base','value',value);
        evalin('base',[fieldname '=value;'])
        
    end
   
    Parallel=myinputs.Parallel;
catch err
display(err)    
end

%%

if nargin<4,
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

MhDirectoryName=myinputs.MhDirectoryName;
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
if ~exist('MhDirectoryName'),
    MhDirectoryName = CheckPath('metropolis',M_.dname);
end

ALPHA = 0.2;                                % increase too much with the number of simulations.
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
    if exist('OCTAVE_VERSION'),
        if (whoiam==0),
            printf('    Parameter %d...  ',j);
        end
    else
        fprintf('    Parameter %d...  ',j);
    end
    for b = 1:nblck
        startline = 0;
        for n = 1:NumberOfMcFilesPerBlock
            %load([MhDirectoryName '/' mcfiles(n,1,b).name],'x2');
            load([MhDirectoryName '/' M_.fname '_mh',int2str(n),'_blck' int2str(b) ...
                  '.mat'],'x2');
            nx2 = size(x2,1);
            tmp((b-1)*NumberOfDraws+startline+(1:nx2),1) = x2(:,j);
            %      clear x2;
            startline = startline + nx2;
        end
% $$$     %load([MhDirectoryName '/' mcfiles(NumberOfMcFilesPerBlock,1,b).name],'x2');
% $$$     load([MhDirectoryName '/' M_.fname '_mh',int2str(NumberOfMcFilesPerBlock),'_blck' int2str(b) '.mat'],'x2');
% $$$     tmp((b-1)*NumberOfDraws+startline+1:(b-1)*NumberOfDraws+MAX_nruns*(LastFileNumber-1)+LastLineNumber,1) = x2(:,j);
% $$$     clear x2;
% $$$     startline = startline + LastLineNumber;
    end
    tmp(:,2) = kron(transpose(1:nblck),ones(NumberOfDraws,1));
    tmp(:,3) = kron(ones(nblck,1),time');
    tmp = sortrows(tmp,1);
    ligne   = 0;
    for iter  = Origin:StepSize:NumberOfDraws
        ligne = ligne+1;
        linea = ceil(mh_drop*iter);
        n     = iter-linea+1;
        cinf  = round(n*ALPHA/2);
        csup  = round(n*(1-ALPHA/2));
        CINF  = round(nblck*n*ALPHA/2);
        CSUP  = round(nblck*n*(1-ALPHA/2));
        temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
        UDIAG(ligne,1,j-fpar+1) = temp(CSUP,1)-temp(CINF,1);
        moyenne = mean(temp(:,1));%% Pooled mean.
        UDIAG(ligne,3,j-fpar+1) = sum((temp(:,1)-moyenne).^2)/(nblck*n-1);
        UDIAG(ligne,5,j-fpar+1) = sum(abs(temp(:,1)-moyenne).^3)/(nblck*n-1);
        for i=1:nblck
            pmet = temp(find(temp(:,2)==i));
            UDIAG(ligne,2,j-fpar+1) = UDIAG(ligne,2,j-fpar+1) + pmet(csup,1)-pmet(cinf,1);
            moyenne = mean(pmet,1); %% Within mean.
            UDIAG(ligne,4,j-fpar+1) = UDIAG(ligne,4,j-fpar+1) + sum((pmet(:,1)-moyenne).^2)/(n-1);
            UDIAG(ligne,6,j-fpar+1) = UDIAG(ligne,6,j-fpar+1) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
        end
    end
    if exist('OCTAVE_VERSION'),
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