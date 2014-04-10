function disp_moments(y,var_list)
% Displays moments of simulated variables

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

global M_ options_ oo_

warning_old_state = warning;
warning off

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
end

nvar = size(var_list,1);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
    if isempty(i_tmp)
        error (['One of the variable specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

y = y(ivar,options_.drop+1:end)';

m = mean(y);

if options_.hp_filter
    [hptrend,y] = sample_hp_filter(y,options_.hp_filter);
else
    y = bsxfun(@minus, y, m);
end

s2 = mean(y.*y);
s = sqrt(s2);
oo_.mean = transpose(m);
oo_.var = y'*y/size(y,1);

labels = deblank(M_.endo_names(ivar,:));

if options_.nomoments == 0
    z = [ m' s' s2' (mean(y.^3)./s2.^1.5)' (mean(y.^4)./(s2.*s2)-3)' ];    
    title='MOMENTS OF SIMULATED VARIABLES';
    if options_.hp_filter
        title = [title ' (HP filter, lambda = ' ...
                 num2str(options_.hp_filter) ')'];
    end
    headers=char('VARIABLE','MEAN','STD. DEV.','VARIANCE','SKEWNESS', ...
                 'KURTOSIS');
    dyntable(title,headers,labels,z,size(labels,2)+2,16,6);
end

if options_.nocorr == 0
    corr = (y'*y/size(y,1))./(s'*s);
    if options_.noprint == 0
        title = 'CORRELATION OF SIMULATED VARIABLES';
        if options_.hp_filter
            title = [title ' (HP filter, lambda = ' ...
                     num2str(options_.hp_filter) ')'];
        end
        headers = char('VARIABLE',M_.endo_names(ivar,:));
        dyntable(title,headers,labels,corr,size(labels,2)+2,8,4);
    end
end

if options_.noprint == 0 && length(options_.conditional_variance_decomposition)
   fprintf('\nSTOCH_SIMUL: conditional_variance_decomposition requires theoretical moments, i.e. periods=0.\n') 
end

ar = options_.ar;
if ar > 0
    autocorr = [];
    for i=1:ar
        oo_.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));
        autocorr = [ autocorr diag(oo_.autocorr{i}) ];
    end
    if options_.noprint == 0
        title = 'AUTOCORRELATION OF SIMULATED VARIABLES';
        if options_.hp_filter
            title = [title ' (HP filter, lambda = ' ...
                     num2str(options_.hp_filter) ')'];
        end
        headers = char('VARIABLE',int2str([1:ar]'));
        dyntable(title,headers,labels,autocorr,size(labels,2)+2,8,4);
    end
end

warning(warning_old_state);
