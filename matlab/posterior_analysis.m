function oo_ = posterior_analysis(type,arg1,arg2,arg3,options_,M_,oo_)
% Copyright (C) 2008-2017 Dynare Team
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

info = check_posterior_analysis_data(type,M_);
SampleSize = options_.sub_draws;
switch info
  case 0
    disp('check_posterior_analysis_data:: Can''t find any mcmc file!')
    error('Check the options of the estimation command...')
  case {1,2}
    MaxMegaBytes = options_.MaximumNumberOfMegaBytes;
    drsize = size_of_the_reduced_form_model(oo_.dr);
    if drsize*SampleSize>MaxMegaBytes
        drsize=0;
    end
    SampleAddress = selec_posterior_draws(SampleSize,drsize);
    oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_);
  case {4,5}
    oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_);
  case 6
    [ivar,vartan] = get_variables_list(options_,M_);
    nvar = length(ivar);
    oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_,nvar,vartan);
  otherwise
    error(['posterior_analysis:: Check_posterior_analysis_data gave a meaningless output!'])
end



function oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_,nvar,vartan)
narg1 = 8;
narg2 = 10;
if ~(nargin==narg1 || nargin==narg2)
    error('posterior_analysis:: Call to function job is buggy!')
end
switch type
  case 'variance'
    if nargin==narg1
        [nvar,vartan,NumberOfFiles] = ...
            dsge_simulated_theoretical_covariance(SampleSize,M_,options_,oo_,'posterior');
    end
    oo_ = covariance_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                 vartan,nvar,arg1,arg2,options_.mh_conf_sig,oo_,options_);
  case 'decomposition'
    if nargin==narg1
        [nvar,vartan,NumberOfFiles] = ...
            dsge_simulated_theoretical_variance_decomposition(SampleSize,M_,options_,oo_,'posterior');
    end
    oo_ = variance_decomposition_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                             M_.exo_names,arg2,vartan,arg1,options_.mh_conf_sig,oo_,options_);
    if ~all(M_.H==0)
        if strmatch(arg1,options_.varobs,'exact')
            [observable_name_requested_vars,index_subset,index_observables]=intersect(vartan,options_.varobs,'stable');            
            oo_ = variance_decomposition_ME_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                             M_.exo_names,arg2,observable_name_requested_vars,arg1,options_.mh_conf_sig,oo_,options_);
        end
    end        
  case 'correlation'
    if nargin==narg1
        [nvar,vartan,NumberOfFiles] = ...
            dsge_simulated_theoretical_correlation(SampleSize,arg3,M_,options_,oo_,'posterior');
    end
    oo_ = correlation_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                  vartan,nvar,arg1,arg2,arg3,options_.mh_conf_sig,oo_,M_,options_);
  case 'conditional decomposition'
    if nargin==narg1
        [nvar,vartan,NumberOfFiles] = ...
            dsge_simulated_theoretical_conditional_variance_decomposition(SampleSize,arg3,M_,options_,oo_,'posterior');
    end
    oo_ = conditional_variance_decomposition_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                                      arg3,M_.exo_names,arg2,vartan,arg1,options_.mh_conf_sig,oo_,options_);
    if ~all(M_.H==0)
        if strmatch(vartan(arg1,:),options_.varobs,'exact')
            oo_ = conditional_variance_decomposition_ME_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                arg3,M_.exo_names,arg2,vartan,arg1,options_.mh_conf_sig,oo_,options_);
        end
    end        
  otherwise
    disp('Not yet implemented')
end