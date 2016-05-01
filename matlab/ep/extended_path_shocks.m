function [shocks, exp_shocks, spfm_exo_simul, innovations, DynareResults] = ...
    extended_path_shocks(innovations, ep, sample_size, DynareModel, DynareResults); 

% Copyright (C) 2016 Dynare Team
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
    
% Simulate shocks.
det_exp_shocks = DynareModel.det_exp_shocks;
det_unexp_shocks = DynareModel.det_unexp_shocks;
if isempty(det_exp_shocks) && isempty(det_unexp_shocks)
    switch ep.innovation_distribution
      case 'gaussian'
        shocks = transpose(transpose(innovations.covariance_matrix_upper_cholesky)*randn(innovations.effective_number_of_shocks,sample_size));
        shocks(:,innovations.positive_var_indx) = shocks;
      otherwise
        error(['extended_path:: ' ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end
else
    shocks = zeros(sample_size,DynareModel.exo_nbr);
    exp_shocks = shocks;
    for i = 1:length(det_exp_shocks)
        k = det_exp_shocks(i).periods;
        ivar = det_exp_shocks(i).exo_id;
        exp_shocks(k,ivar) = det_exp_shocks(i).value;
    end

    for i = 1:length(det_unexp_shocks)
        k = det_unexp_shocks(i).periods;
        ivar = det_unexp_shocks(i).exo_id;
        shocks(k,ivar) = det_unexp_shocks(i).value;
    end

    innovations.positive_var_indx = find(sum(abs(shocks)>0));
    innovations.exp_positive_var_indx = find(sum(abs(exp_shocks)>0));
end

% Copy the shocks in exo_simul
DynareResults.exo_simul = shocks;
DynareResults.exp_exo_simul = exp_shocks;
spfm_exo_simul = repmat(DynareResults.exo_steady_state',ep.periods+2,1);