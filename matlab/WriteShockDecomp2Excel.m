function WriteShockDecomp2Excel(z,shock_names,endo_names,i_var,initial_date,DynareModel,DynareOptions,opts_decomp)
%function WriteShockDecomp2Excel(z,shock_names,endo_names,i_var,initial_date,DynareModel,DynareOptions)
% Saves the results from the shock_decomposition command to xls
% 
% Inputs
%   z               [n_var*(nshock+2)*nperiods]     shock decomposition array, see shock_decomposition.m for details
%   shock_names     [endo_nbr*string length]        shock names from M_.exo_names
%   endo_names      [exo_nbr*string length]         variable names from M_.endo_names
%   i_var           [n_var*1]                       vector indices of requested variables in M_.endo_names and z
%   initial_date    [dseries object]                first period of decomposition to plot
%   DynareModel     [structure]                     Dynare model structure
%   DynareOptions   [structure]                     Dynare options structure

% Copyright (C) 2016-2017 Dynare Team
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

SteadyState=[];
fig_mode='';
fig_mode1='';
fig_names='';
screen_shocks=0;
use_shock_groups = DynareOptions.use_shock_groups;
if use_shock_groups
    shock_groups = DynareModel.shock_groups.(use_shock_groups);
    shock_ind = fieldnames(shock_groups);
end

% number of components equals number of shocks + 1 (initial conditions)
comp_nbr = size(z,2)-1;

if nargin==8 ,    
    if isfield(opts_decomp,'steady_state')
        SteadyState = opts_decomp.steady_state;
    end
    if isfield(opts_decomp,'fig_mode') && ~isempty(opts_decomp.fig_mode)
        fig_mode = opts_decomp.fig_mode;
        fig_mode1 = ['_' fig_mode];        
        fig_mode = [fig_mode '_'];        
    end
    if isfield(opts_decomp,'screen_shocks')
        if DynareOptions.use_shock_groups,
            screen_shocks=0;
        elseif comp_nbr>18
            screen_shocks = opts_decomp.screen_shocks;
        end
    end
    if isfield(opts_decomp,'fig_names')
        fig_names = opts_decomp.fig_names;
%         fig_names = ['_' fig_names];
        fig_names1 = [fig_names];        
        fig_names = [fig_names '_'];    
    end
    if screen_shocks
        fig_names1 = [fig_names1 '_screen'];        
        fig_names = [fig_names 'screen_'];            
    end
end   


gend = size(z,3);
if isempty(initial_date)
    x = 1:gend;
else
    freq = initial_date.freq;
    initial_period = initial_date.time(1) + (initial_date.time(2)-1)/freq;
    x = initial_period:(1/freq):initial_period+(gend-1)/freq;
end


nvar = length(i_var);

labels = char(char(shock_names),'Initial values');
if ~(screen_shocks && comp_nbr>18),
    screen_shocks=0;
end
comp_nbr0=comp_nbr;
%%plot decomposition
for j=1:nvar
    d0={};
    z1 = squeeze(z(i_var(j),:,:));
    if screen_shocks,
        [junk, isort] = sort(mean(abs(z1(1:end-2,:)')), 'descend');
        labels = char(char(shock_names(isort(1:16),:)),'Others', 'Initial values');
        zres = sum(z1(isort(17:end),:),1);
        z1 = [z1(isort(1:16),:); zres; z1(comp_nbr0:end,:)];
        comp_nbr=18;
    end
    
    d0(1,:)=[{'Decomposition'} cellstr(labels(1:comp_nbr,:))' {'Smoot Var'}];
    d0=[d0; num2cell([x' z1'])];
    LastRow=size(d0,1);
    if use_shock_groups
        d0(LastRow+2,1)={'Legend.'};
        d0(LastRow+2,2)={'Shocks include:'};
        d0(LastRow+3:LastRow+3+comp_nbr-1,1)=cellstr(labels(1:comp_nbr,:));
        for ic=1:comp_nbr,
            group_members = shock_groups.(shock_ind{ic}).shocks;
            d0(LastRow+2+ic,2:1+length(group_members))=group_members;
        end
    end
    
    warning off
    if ~ismac
        [STATUS,MESSAGE] = xlswrite([DynareModel.fname,'_shock_decomposition',fig_mode,fig_names1],d0,deblank(endo_names(i_var(j),:)));
    else
        [STATUS] = xlwrite([DynareModel.fname,'_shock_decomposition',fig_mode,fig_names1],d0,deblank(endo_names(i_var(j),:)));
    end
    warning on

    clear d0
    
end

