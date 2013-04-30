function [omega,f] = UnivariateSpectralDensity(dr,var_list)
% This function computes the theoretical spectral density of each
% endogenous variable declared in var_list. Results are stored in 
% oo_ and may be plotted. Plots are saved into the graphs-folder.
% 
% Adapted from th_autocovariances.m.  

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

global options_ oo_ M_


omega = []; f = [];

if options_.order > 1
    disp('UnivariateSpectralDensity :: I Cannot compute the theoretical spectral density') 
    disp('with a second order approximation of the DSGE model!')
    disp('Please set order = 1. I abort')
    return
end
pltinfo  = options_.SpectralDensity.plot; 
cutoff   = options_.SpectralDensity.cutoff; 
sdl      = options_.SpectralDensity.sdl; 
omega    = (0:sdl:pi)';
GridSize = length(omega);
exo_names_orig_ord  = M_.exo_names_orig_ord;
if exist('OCTAVE_VERSION')
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
f = zeros(nvar,GridSize);
ghx = dr.ghx;
ghu = dr.ghu;
nspred = M_.nspred;
nstatic = M_.nstatic;
kstate = dr.kstate;
order = dr.order_var;
iv(order) = [1:length(order)];
nx = size(ghx,2);
ikx = [nstatic+1:nstatic+nspred];
A = zeros(nx,nx);
k0 = kstate(find(kstate(:,2) <= M_.maximum_lag+1),:);
i0 = find(k0(:,2) == M_.maximum_lag+1);
i00 = i0;
n0 = length(i0);
A(i0,:) = ghx(ikx,:);
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
Gamma = zeros(nvar,cutoff+1);
[A,B] = kalman_transition_matrix(dr,ikx',1:nx,M_.exo_nbr);
[vx, u] =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold);
iky = iv(ivar);
if ~isempty(u)
    iky = iky(find(any(abs(ghx(iky,:)*u) < options_.Schur_vec_tol,2)));
    ivar = dr.order_var(iky);
end
aa = ghx(iky,:);
bb = ghu(iky,:);

if options_.hp_filter == 0
    tmp = aa*vx*aa'+ bb*M_.Sigma_e*bb';
    k = find(abs(tmp) < 1e-12);
    tmp(k) = 0;
    Gamma(:,1) = diag(tmp);
    vxy = (A*vx*aa'+ghu1*M_.Sigma_e*bb');
    tmp = aa*vxy;
    k = find(abs(tmp) < 1e-12);
    tmp(k) = 0;
    Gamma(:,2) = diag(tmp);
    for i=2:cutoff
        vxy = A*vxy;
        tmp = aa*vxy;
        k = find(abs(tmp) < 1e-12);
        tmp(k) = 0;
        Gamma(:,i+1) = diag(tmp);
    end
else
    iky = iv(ivar);  
    aa = ghx(iky,:);
    bb = ghu(iky,:);
    lambda = options_.hp_filter;
    ngrid = options_.hp_ngrid;
    freqs = 0 : ((2*pi)/ngrid) : (2*pi*(1 - .5/ngrid)); 
    tpos  = exp( sqrt(-1)*freqs);
    tneg  =  exp(-sqrt(-1)*freqs);
    hp1 = 4*lambda*(1 - cos(freqs)).^2 ./ (1 + 4*lambda*(1 - cos(freqs)).^2);
    mathp_col = [];
    IA = eye(size(A,1));
    IE = eye(M_.exo_nbr);
    for ig = 1:ngrid
        f_omega  =(1/(2*pi))*( [inv(IA-A*tneg(ig))*ghu1;IE]...
                               *M_.Sigma_e*[ghu1'*inv(IA-A'*tpos(ig)) IE]); % state variables
        g_omega = [aa*tneg(ig) bb]*f_omega*[aa'*tpos(ig); bb']; % selected variables
        f_hp = hp1(ig)^2*g_omega; % spectral density of selected filtered series
        mathp_col = [mathp_col ; (f_hp(:))'];    % store as matrix row
                                                 % for ifft
    end;
    imathp_col = real(ifft(mathp_col))*(2*pi);
    tmp = reshape(imathp_col(1,:),nvar,nvar);
    k = find(abs(tmp)<1e-12);
    tmp(k) = 0;
    Gamma(:,1) = diag(tmp);
    sy = sqrt(Gamma(:,1));
    sy = sy *sy';
    for i=1:cutoff-1
        tmp = reshape(imathp_col(i+1,:),nvar,nvar)./sy;
        k = find(abs(tmp) < 1e-12);
        tmp(k) = 0;
        Gamma(:,i+1) = diag(tmp);
    end
end

H = 1:cutoff;
for i=1:nvar
    f(i,:) = Gamma(i,1)/(2*pi) + Gamma(i,H+1)*cos(H'*omega')/pi;
end  

if exist('OCTAVE_VERSION')
    warning('on', 'Octave:divide-by-zero')
else
    warning on MATLAB:dividebyzero
end

if pltinfo
    if ~exist(M_.fname, 'dir')
        mkdir('.',M_.fname);
    end
    if ~exist([M_.fname '/graphs'])
        mkdir(M_.fname,'graphs');
    end

    for i= 1:nvar
        hh = dyn_figure(options_,'Name',['Spectral Density of ' deblank(M_.endo_names(ivar(i),:)) '.']);
        plot(omega,f(i,:),'-k','linewidth',2)
        xlabel('0 \leq \omega \leq \pi')
        ylabel('f(\omega)')
        box on
        axis tight        
        dyn_saveas(hh,[M_.fname ,filesep,'graphs', filesep, 'SpectralDensity_' deblank(M_.endo_names(ivar(i),:))],options_)
    end
end
