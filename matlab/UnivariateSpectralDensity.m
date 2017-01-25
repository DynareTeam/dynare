function [oo_] = UnivariateSpectralDensity(M_,oo_,options_,var_list)
% This function computes the theoretical spectral density of each
% endogenous variable declared in var_list. Results are stored in 
% oo_.SpectralDensity and may be plotted. Plots are saved into the 
% graphs-folder. 
% 
% INPUTS
%   M_                  [structure]    Dynare's model structure
%   oo_                 [structure]    Dynare's results structure
%   options_            [structure]    Dynare's options structure
%   var_list            [integer]      Vector of indices for a subset of variables.
%    
% OUTPUTS
%   oo_                 [structure]    Dynare's results structure,
%                                       containing the subfield
%                                       SpectralDensity with fields freqs
%                                       and density, which are of size nvar*ngrid.
%

% Adapted from th_autocovariances.m.  

% Copyright (C) 2006-2015 Dynare Team
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


if options_.order > 1
    disp('UnivariateSpectralDensity :: I Cannot compute the theoretical spectral density') 
    disp('with a second order approximation of the DSGE model!')
    disp('Please set order = 1. I abort')
    return
end

if isoctave
    warning('off', 'Octave:divide-by-zero')
else
    warning off MATLAB:dividebyzero
end
if nargin<2
    var_list = [];
end
if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
end
nvar = size(var_list,1);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
    if isempty(i_tmp)
        error (['One of the variables specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

ghx = oo_.dr.ghx;
ghu = oo_.dr.ghu;
nspred = M_.nspred;
nstatic = M_.nstatic;
kstate = oo_.dr.kstate;
order = oo_.dr.order_var;
iv(order) = [1:length(order)];
nx = size(ghx,2);
ikx = [nstatic+1:nstatic+nspred];
k0 = kstate(find(kstate(:,2) <= M_.maximum_lag+1),:);
i0 = find(k0(:,2) == M_.maximum_lag+1);
i00 = i0;
AS = ghx(:,i0);
ghu1 = zeros(nx,M_.exo_nbr);
ghu1(i0,:) = ghu(ikx,:);
for i=M_.maximum_lag:-1:2
    i1 = find(k0(:,2) == i);
    n1 = size(i1,1);
    j1 = zeros(n1,1);
    j2 = j1;
    for k1 = 1:n1
        j1(k1) = find(k0(i00,1)==k0(i1(k1),1));
        j2(k1) = find(k0(i0,1)==k0(i1(k1),1));
    end
    AS(:,j1) = AS(:,j1)+ghx(:,i1);
    i0 = i1;
end

[A,B] = kalman_transition_matrix(oo_.dr,ikx',1:nx,M_.exo_nbr);
[vx, u] =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
iky = iv(ivar);
if ~isempty(u)
    iky = iky(find(any(abs(ghx(iky,:)*u) < options_.Schur_vec_tol,2)));
    ivar = oo_.dr.order_var(iky);
end

iky = iv(ivar);  
aa = ghx(iky,:);
bb = ghu(iky,:);
ngrid = options_.hp_ngrid; %number of grid points
freqs = (0 : pi/(ngrid-1):pi)'; % grid on which to compute
tpos  = exp( sqrt(-1)*freqs); %positive frequencies
tneg  = exp(-sqrt(-1)*freqs); %negative frequencies
if options_.one_sided_hp_filter
    error('UnivariateSpectralDensity:: spectral density estimate not available with one-sided HP filter')
elseif options_.hp_filter == 0 && ~options_.bandpass.indicator %do not filter
   filter_gain=ones(ngrid,1);
elseif ~(options_.hp_filter == 0 && ~options_.bandpass.indicator) && options_.bandpass.indicator %filter with bandpass
    filter_gain = zeros(1,ngrid);
    lowest_periodicity=options_.bandpass.passband(2);
    highest_periodicity=options_.bandpass.passband(1);
    highest_periodicity=max(2,highest_periodicity); % restrict to upper bound of pi
    filter_gain(freqs>=2*pi/lowest_periodicity & freqs<=2*pi/highest_periodicity)=1;
    filter_gain(freqs<=-2*pi/lowest_periodicity+2*pi & freqs>=-2*pi/highest_periodicity+2*pi)=1;
elseif ~(options_.hp_filter == 0 && ~options_.bandpass.indicator) && ~options_.bandpass.indicator %filter with HP-filter
    lambda = options_.hp_filter;
    filter_gain = 4*lambda*(1 - cos(freqs)).^2 ./ (1 + 4*lambda*(1 - cos(freqs)).^2);    
end

mathp_col = NaN(ngrid,length(ivar)^2);
IA = eye(size(A,1));
IE = eye(M_.exo_nbr);
for ig = 1:ngrid
    f_omega  =(1/(2*pi))*( [(IA-A*tneg(ig))\ghu1;IE]...
                           *M_.Sigma_e*[ghu1'/(IA-A'*tpos(ig)) IE]); % state variables
    g_omega = [aa*tneg(ig) bb]*f_omega*[aa'*tpos(ig); bb']; % selected variables
    f_hp = filter_gain(ig)^2*g_omega; % spectral density of selected filtered series
    mathp_col(ig,:) = (f_hp(:))';    % store as matrix row
end;

f = zeros(nvar,ngrid);
for i=1:nvar
    f(i,:) = real(mathp_col(:,(i-1)*nvar+i)); %read out spectral density
end  

oo_.SpectralDensity.freqs=freqs;
oo_.SpectralDensity.density=f;

if isoctave
    warning('on', 'Octave:divide-by-zero')
else
    warning on MATLAB:dividebyzero
end

if options_.nograph == 0
    if ~exist(M_.fname, 'dir')
        mkdir('.',M_.fname);
    end
    if ~exist([M_.fname '/graphs'],'dir')
        mkdir(M_.fname,'graphs');
    end

    for i= 1:nvar
        hh = dyn_figure(options_,'Name',['Spectral Density of ' deblank(M_.endo_names(ivar(i),:)) '.']);
        plot(freqs,f(i,:),'-k','linewidth',2)
        xlabel('0 \leq \omega \leq \pi')
        ylabel('f(\omega)')
        box on
        axis tight        
        dyn_saveas(hh,[M_.fname ,filesep,'graphs', filesep, 'SpectralDensity_' deblank(M_.endo_names(ivar(i),:))],options_)
    end
end
