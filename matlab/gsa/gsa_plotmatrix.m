function gsa_plotmatrix(type,varargin)
% function gsa_plotmatrix(type,varargin)
% extended version of the standard MATLAB plotmatrix 
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2011-2012 European Commission
% Copyright (C) 2011-2012 Dynare Team
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

global bayestopt_ options_ M_

RootDirectoryName = CheckPath('gsa',M_.dname);

if options_.opt_gsa.pprior
    load([ RootDirectoryName filesep  M_.fname '_prior.mat'],'lpmat0','lpmat','istable','iunstable','iindeterm','iwrong')
else
    load([ RootDirectoryName filesep  M_.fname '_mc.mat'],'lpmat0','lpmat','istable','iunstable','iindeterm','iwrong')
    eval(['load ' options_.mode_file ' xparam1;']');
end

iexplosive = iunstable(~ismember(iunstable,[iindeterm;iwrong]));

switch type
    case 'all'
        x=[lpmat0 lpmat];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'stable'
        x=[lpmat0(istable,:) lpmat(istable,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'nosolution'
        x=[lpmat0(iunstable,:) lpmat(iunstable,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'unstable'
        x=[lpmat0(iexplosive,:) lpmat(iexplosive,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'indeterm'
        x=[lpmat0(iindeterm,:) lpmat(iindeterm,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'wrong'
        x=[lpmat0(iwrong,:) lpmat(iwrong,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
        
end

if isempty(x),
    disp('Empty parameter set!')
    return
end

for j=1:length(varargin),
    jcol(j)=strmatch(varargin{j},bayestopt_.name,'exact');
end

[H,AX,BigA,P,PAx]=plotmatrix(x(:,jcol));

for j=1:length(varargin),
    %      axes(AX(1,j)), title(varargin{j})
    %      axes(AX(j,1)), ylabel(varargin{j})
    %      set(AX(1,j),'title',varargin{j}),
    set(get(AX(j,1),'ylabel'),'string',varargin{j})
    set(get(AX(end,j),'xlabel'),'string',varargin{j})
end

if options_.opt_gsa.pprior==0,
    xparam1=xparam1(jcol);
    for j=1:length(varargin),
        for i=1:j-1,
            axes(AX(j,i)),
            hold on, plot(xparam1(i),xparam1(j),'*r')
        end
        for i=j+1:length(varargin),
            axes(AX(j,i)),
            hold on, plot(xparam1(i),xparam1(j),'*r')
        end
    end
end
