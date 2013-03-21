function [nam,texnam] = get_the_name(k,TeX,M_,estim_params_,options_)

%@info:
%! @deftypefn {Function File} {[@var{nam},@var{texnam}] =} get_the_name (@var{k},@var{TeX},@var{M_},@var{estim_params_},@var{options_})
%! @anchor{get_the_name}
%! @sp 1
%! Returns the name of the estimated parameter number @var{k}, following the internal ordering of the estimated parameters.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item k
%! Integer scalar, parameter number.
%! @item TeX
%! Integer scalar, if @var{TeX}==0 then @var{texnam} is not returned (empty matrix).
%! @item M_
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item estim_params_
%! Matlab's structure describing the estimated parameters (initialized by @code{dynare}).
%! @item options_
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item nam
%! String, internal name of the variable
%! @item texnam
%! String, TeX name of the same variable (if defined in the mod file).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{get_prior_info}, @ref{McMCDiagnostics}, @ref{mode_check}, @ref{PlotPosteriorDistributions}, @ref{plot_priors}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! None.
%! @end deftypefn
%@eod:

% Copyright (C) 2004-2013 Dynare Team
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

nam = [];
texnam = [];

nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;

if k <= nvx
    vname = deblank(M_.exo_names(estim_params_.var_exo(k,1),:));
    nam = ['SE_',vname];
    if TeX
        tname  = deblank(M_.exo_names_tex(estim_params_.var_exo(k,1),:));
        texnam = ['$ SE_{' tname '} $'];
    end
elseif  k <= (nvx+nvn)
    vname = deblank(options_.varobs(estim_params_.nvn_observable_correspondence(k-estim_params_.nvx,1),:));
    nam=['SE_EOBS_',vname];
    if TeX
        tname  = deblank(M_.endo_names_tex(estim_params_.var_endo(k-estim_params_.nvx,1),:));
        texnam = ['$ EOBS SE_{' tname '} $'];
    end
elseif  k <= (nvx+nvn+ncx)
    jj = k - (nvx+nvn);
    k1 = estim_params_.corrx(jj,1);
    k2 = estim_params_.corrx(jj,2);
    vname = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
    nam=['CC_',vname];
    if TeX
        tname  = [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))];
        texnam = ['$ CC_{' tname '} $'];
    end
elseif  k <= (nvx+nvn+ncx+ncn)
    jj = k - (nvx+nvn+ncx);
    k1 = estim_params_.corrn(jj,1);
    k2 = estim_params_.corrn(jj,2);
    vname = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
    nam=['CC_EOBS_' vname];
    if TeX
        tname  = [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))];
        texnam =['$ EOBS CC_{' tname '} $'];
    end
else
    jj = k - (nvx+nvn+ncx+ncn);
    jj1 = estim_params_.param_vals(jj,1);
    nam = deblank(M_.param_names(jj1,:));
    if TeX
        texnam = ['$ '  deblank(M_.param_names_tex(jj1,:))  ' $'];
    end
end