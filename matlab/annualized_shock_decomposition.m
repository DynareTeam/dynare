function [z, endo_names, endo_names_tex, steady_state, i_var, oo_] = annualized_shock_decomposition(oo_, M_, options_, i_var, t0, t1, realtime_, vintage_, steady_state, q2a, cumfix)
% function oo_ = annualized_shock_decomposition(oo_,t0,options_.nobs);
% Computes annualized shocks contribution to a simulated trajectory. The fields set are
% oo_.annualized_shock_decomposition, oo_.annualized_realtime_shock_decomposition, 
% oo_.annualized_realtime_conditional_shock_decomposition and oo_.annualized_realtime_forecast_shock_decomposition. 
% Subfields are arrays n_var by nshock+2 by nperiods. The
% first nshock columns store the respective shock contributions, column n+1
% stores the role of the initial conditions, while column n+2 stores the
% value of the smoothed variables.  Both the variables and shocks are stored 
% in the order of endo_names and M_.exo_names, respectively.
%
% INPUTS
%    oo_:          [structure] Storage of results
%    M_:           [structure] Storage of model
%    opts:         [structure] options for shock decomp 
%    i_var:        [array] index of vars
%    t0:           [integer]  first period
%    t1:           [integer] last period
%    realtime_:    [integer]   
%    vintage_:     [integer]
%    steady_state: [array] steady state value of quarterly (log-) level vars
%    q2a:          [structure] info on q2a
%
% OUTPUTS
%    z:              [matrix] shock decomp to plot 
%    endo_names:     [char] updated var names
%    endo_names_tex: [char] updated TeX var names
%    steady_state:   [array] updated stady state of vars
%    i_var:          [integer array] updated var indices to plot
%    oo_:            [structure]  Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017 Dynare Team
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

opts = options_.shock_decomp;
nvar = length(i_var);
GYTREND0 = q2a.GYTREND0;
var_type = q2a.type;
islog    = q2a.islog;
aux      = q2a.aux;
aux0 = aux;
cumfix      = q2a.cumfix;
% usual shock decomp 
if isstruct(oo_)
%     z = oo_.shock_decomposition;
        myopts=options_;
        myopts.shock_decomp.type='qoq';
        myopts.shock_decomp.realtime=0;
        [z, junk] = plot_shock_decomposition(M_,oo_,myopts,[]);
else
    z = oo_;
end
z = z(i_var,:,:);
mytype=var_type;
if isfield(q2a,'name')
    mytxt = q2a.name;
    mytex = q2a.name;
    if isfield(q2a,'tex_name')
        mytex = q2a.tex_name;
    end
    if mytype==2,
        gtxt = ['PHI' mytxt]; % inflation rate
        gtex = ['{\pi(' mytex ')}'];
    elseif mytype
        gtxt = ['G' mytxt]; % inflation rate
        gtex = ['{g(' mytex ')}'];
    end
    if isfield(q2a,'gname')
        gtxt = q2a.gname;
    end
    if isfield(q2a,'tex_gname')
        gtex = q2a.tex_gname;
    end
    mytype=0;
end
if isstruct(aux)
    if ischar(aux.y)
        myopts=options_;
        myopts.shock_decomp.type='qoq';
        myopts.shock_decomp.realtime=0;
        [y_aux, steady_state_aux] = plot_shock_decomposition(M_,oo_,myopts,aux.y);
        aux.y=y_aux;
        aux.yss=steady_state_aux;
    end
    yaux=aux.y;
end
if mytype==2,
    gtxt = 'PHI'; % inflation rate
    gtex = '\pi';
elseif mytype
    gtxt = 'G'; % growth rate
    gtex = 'g';
end
steady_state=steady_state(i_var);
% endo_names = M_.endo_names(i_var,:);
% endo_names_tex = M_.endo_names_tex(i_var,:);
nterms = size(z,2);
nfrcst = opts.forecast/4;

for j=1:nvar
    if j>1,
        endo_names = char(endo_names,[deblank(M_.endo_names(i_var(j),:)) '_A']);
        endo_names_tex = char(endo_names_tex,['{' deblank(M_.endo_names_tex(i_var(j),:)) '}^A']);
        gendo_names = char(gendo_names,[gtxt endo_names(j,:)]);
        gendo_names_tex = char(gendo_names_tex,[gtex '(' deblank(endo_names_tex(j,:)) ')']);
    else
        if nvar==1 && ~mytype
            endo_names = mytxt;
            endo_names_tex = mytex;
            gendo_names = gtxt;
            gendo_names_tex = gtex;
        else
            endo_names = [deblank(M_.endo_names(i_var(j),:)) '_A'];
            endo_names_tex = ['{' deblank(M_.endo_names_tex(i_var(j),:)) '}^A'];
            gendo_names = [gtxt endo_names(j,:)];
            gendo_names_tex = [gtex '(' deblank(endo_names_tex(j,:)) ')'];
        end
    end
    for k =1:nterms,
        if isstruct(aux),
            aux.y = squeeze(yaux(j,k,min((t0-3):-4:1):end));
        end
        [za(j,k,:), steady_state_a(j,1), gza(j,k,:), steady_state_ga(j,1)] = ...
            quarterly2annual(squeeze(z(j,k,min((t0-3):-4:1):end)),steady_state(j),GYTREND0,var_type,islog,aux);
    end
    ztmp=squeeze(za(j,:,:));
    if cumfix==0,
        zscale = sum(ztmp(1:end-1,:))./ztmp(end,:);
        ztmp(1:end-1,:) = ztmp(1:end-1,:)./repmat(zscale,[nterms-1,1]);
    else
        zres = ztmp(end,:)-sum(ztmp(1:end-1,:));
        ztmp(end-1,:) = ztmp(end-1,:) + zres;
    end
    gztmp=squeeze(gza(j,:,:));
    if cumfix==0,
        gscale = sum(gztmp(1:end-1,:))./ gztmp(end,:);
        gztmp(1:end-1,:) = gztmp(1:end-1,:)./repmat(gscale,[nterms-1,1]);
    else
        gres = gztmp(end,:) - sum(gztmp(1:end-1,:));
        gztmp(end-1,:) = gztmp(end-1,:)+gres;
    end
    za(j,:,:) = ztmp;
    gza(j,:,:) = gztmp;
end

if q2a.plot ==1,
    z=gza;
    endo_names = gendo_names;
    endo_names_tex = gendo_names_tex;
elseif q2a.plot == 2
    z=za;
else
    z=cat(1,za,gza);
    endo_names = char(endo_names,gendo_names);
    endo_names_tex = char(endo_names_tex,gendo_names_tex);
end
% if isstruct(oo_)
%     oo_.annualized_shock_decomposition=z;
% end

% realtime
if realtime_ && isstruct(oo_) && isfield(oo_, 'realtime_shock_decomposition'),
init=1;
for i=t0:4:t1,
    yr=floor(i/4);
    za=[];
    gza=[];
        myopts=options_;
        myopts.shock_decomp.type='qoq';
        myopts.shock_decomp.realtime=1;
        myopts.shock_decomp.vintage=i;
        [z, steady_state_aux] = plot_shock_decomposition(M_,oo_,myopts,[]);
        z = z(i_var,:,:);
if isstruct(aux)
    if ischar(aux0.y)
        [y_aux, steady_state_aux] = plot_shock_decomposition(M_,oo_,myopts,aux0.y);
        aux.y=y_aux;
        aux.yss=steady_state_aux;
    end
    yaux=aux.y;
end
        nterms = size(z,2);
  
%     z = oo_.realtime_shock_decomposition.(['time_' int2str(i)]);
%     z = z(i_var,:,:);
           
    for j=1:nvar
        for k =nterms:-1:1,
%             if k<nterms
%                 ztmp = squeeze(sum(z(j,[1:k-1,k+1:end-1],t0-4:end)));
%             else
                ztmp = squeeze(z(j,k,min((t0-3):-4:1):end));
%             end
            if isstruct(aux),
                aux.y = squeeze(yaux(j,k,min((t0-3):-4:1):end));
            end
            [za(j,k,:), steady_state_a(j,1), gza(j,k,:), steady_state_ga(j,1)] = ...
                quarterly2annual(ztmp,steady_state(j),GYTREND0,var_type,islog,aux);
%             if k<nterms
%                 za(j,k,:) = za(j,end,:) - za(j,k,:);
%                 gza(j,k,:) = gza(j,end,:) - gza(j,k,:);
%             end
            
        end
        
        ztmp=squeeze(za(j,:,:));

        if cumfix==0,
            zscale = sum(ztmp(1:end-1,:))./ztmp(end,:);
            ztmp(1:end-1,:) = ztmp(1:end-1,:)./repmat(zscale,[nterms-1,1]);
        else
            zres = ztmp(end,:)-sum(ztmp(1:end-1,:));
            ztmp(end-1,:) = ztmp(end-1,:) + zres;
        end
        
        gztmp=squeeze(gza(j,:,:));
        if cumfix==0,
            gscale = sum(gztmp(1:end-1,:))./ gztmp(end,:);
            gztmp(1:end-1,:) = gztmp(1:end-1,:)./repmat(gscale,[nterms-1,1]);
        else
            gres = gztmp(end,:) - sum(gztmp(1:end-1,:));
            gztmp(end-1,:) = gztmp(end-1,:)+gres;
        end
        
        za(j,:,:) = ztmp;
        gza(j,:,:) = gztmp;
    end
    
    if q2a.plot ==1,
        z=gza;
    elseif q2a.plot == 2
        z=za;
    else
        z=cat(1,za,gza);
    end
    
    if init==1,
        oo_.annualized_realtime_shock_decomposition.pool = z;
    else
        oo_.annualized_realtime_shock_decomposition.pool(:,:,yr) = z(:,:,end-nfrcst);
    end        
    oo_.annualized_realtime_shock_decomposition.(['yr_' int2str(yr)]) = z;
        
    if opts.forecast
        oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(yr)]) = z(:,:,end-nfrcst:end);
        if init>nfrcst
            oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(yr-nfrcst)]) = ...
                oo_.annualized_realtime_shock_decomposition.pool(:,:,yr-nfrcst:end) - ...
                oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(yr-nfrcst)]);
            % fix others
            oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(yr-nfrcst)])(:,end-1,:) = ...
                oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(yr-nfrcst)])(:,end-1,:) + ...
                oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(yr-nfrcst)])(:,end,:);
            % fix total
            oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(yr-nfrcst)])(:,end,:) = ...
                oo_.annualized_realtime_shock_decomposition.pool(:,end,yr-nfrcst:end);
        end
    end
% ztmp=oo_.realtime_shock_decomposition.pool(:,:,21:29)-oo_.realtime_forecast_shock_decomposition.time_21;



    init=init+1;
end


switch realtime_
    
    case 0
        z = oo_.annualized_shock_decomposition;
    
    case 1 % realtime
        if vintage_
            z = oo_.annualized_realtime_shock_decomposition.(['yr_' int2str(floor(vintage_/4))]);
        else
            z = oo_.annualized_realtime_shock_decomposition.pool;
        end
    
    case 2 % conditional
        if vintage_
            z = oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(floor(vintage_/4))]);
        else
            error();
        end
        
    case 3 % forecast
        if vintage_
            z = oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(floor(vintage_/4))]);
        else
            error()
        end
end
end

if q2a.plot ==0,
    i_var=1:2*nvar;
    steady_state = [steady_state_a;steady_state_ga];
else
    i_var=1:nvar;
    if q2a.plot ==1,
        steady_state = steady_state_ga;
    else
        steady_state = steady_state_a;
    end
end

