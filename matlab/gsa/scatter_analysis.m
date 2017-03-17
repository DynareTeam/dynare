function indmcf = scatter_analysis(lpmat, xdata, options_scatter, DynareOptions)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu
%

% Copyright (C) 2017 European Commission
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

param_names = options_scatter.param_names;

if DynareOptions.TeX
    if ~isfield(options_scatter,'param_names_tex')
        param_names_tex = options_scatter.param_names;
    else
        param_names_tex = options_scatter.param_names_tex;
    end
end
amcf_name = options_scatter.amcf_name;
amcf_title = options_scatter.amcf_title;
title = options_scatter.title;
fname_ = options_scatter.fname_;
xparam1=[];
if isfield(options_scatter,'xparam1'),
    xparam1=options_scatter.xparam1;
end    
OutputDirectoryName = options_scatter.OutputDirectoryName;

if ~DynareOptions.nograph,
    skipline()
    xx=[];
    if ~ isempty(xparam1), xx=xparam1; end
    scatter_plots(lpmat, xdata, param_names, ...
        '.', [fname_,'_',amcf_name], OutputDirectoryName, amcf_title,xx, DynareOptions)
end
