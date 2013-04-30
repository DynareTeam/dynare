function forecast_graphs(var_list)

% Copyright (C) 2008-2012 Dynare Team
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

global M_ oo_ options_

nc = 4;
nr = 3;
exo_nbr = M_.exo_nbr;
endo_names = M_.endo_names;
fname = M_.fname;
% $$$     varobs = options_.varobs;
% $$$     y = oo_.SmoothedVariables;
% $$$     ys = oo_.dr.ys;
% $$$     gend = size(y,2);
yf = oo_.forecast.Mean;
hpdinf = oo_.forecast.HPDinf;
hpdsup = oo_.forecast.HPDsup;
if isempty(var_list)
    varlist = endo_names(1:M_.orig_endo_nbr,:);
end
i_var = [];
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t and endogenous variable'])
    end
    i_var = [i_var; tmp];
end
nvar = length(i_var);

% $$$     % build trend for smoothed variables if necessary
% $$$     trend = zeros(size(varobs,1),10);
% $$$     if isfield(oo_.Smoother,'TrendCoeffs')
% $$$         trend_coeffs = oo_.Smoother.TrendCoeffs;
% $$$         trend = trend_coeffs*(gend-9:gend);
% $$$     end

% create subdirectory <fname>/graphs if id doesn't exist
if ~exist(fname, 'dir')
    mkdir('.',fname);
end
if ~exist([fname '/graphs'])
    mkdir(fname,'graphs');
end

m = 1;
n_fig = 1;
hh=dyn_figure(options_,'Name','Forecasts (I)');
for j= 1:nvar
    if m > nc*nr; 
        dyn_saveas(hh,[ fname '/graphs/forcst' int2str(n_fig)],options_);
        
        n_fig =n_fig+1;
        eval(['hh=dyn_figure(options_,''Name'',''Forecasts (' int2str(n_fig) ')'');']);
        m = 1;
    end
    subplot(nr,nc,m);
    vn = deblank(endo_names(i_var(j),:));
    obs = 0;
% $$$         k = strmatch(vn,varobs,'exact'); 
% $$$   if ~isempty(k)
% $$$       yy = y.(vn)(end-9:end) + repmat(ys(i_var(j)),10,1)+trend(k,:)';
% $$$       plot(yy);
% $$$       hold on
% $$$       obs = 10;
% $$$   end
    plot([NaN(obs,1); yf.(vn)]);
    hold on
    plot([NaN(obs,1); hpdinf.(vn)]);
    hold on
    plot([NaN(obs,1); hpdsup.(vn)]);
    title(vn,'Interpreter','none');
    xlim([1 obs+length(hpdsup.(vn))])
    hold off
    m = m + 1;
end

if m > 1
    dyn_saveas(hh,[fname '/graphs/forcst' int2str(n_fig)],options_);
end