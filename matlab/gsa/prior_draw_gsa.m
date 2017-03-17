function pdraw = prior_draw_gsa(init,rdraw)
% Draws from the prior distributions for use with Sensitivity Toolbox for DYNARE
% 
% INPUTS
%   o init           [integer]  scalar equal to 1 (first call) or 0.
%   o rdraw          
%    
% OUTPUTS 
%   o pdraw          [double]   draw from the joint prior density. 
%
% ALGORITHM 
%   ...       
%
% SPECIAL REQUIREMENTS
%   MATLAB Statistics Toolbox
%  
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2012-2015 European Commission
% Copyright (C) 2012-2015 Dynare Team
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

global bayestopt_ options_ estim_params_ M_
persistent npar pshape p6 p7 p3 p4 lbcum ubcum
  
if init
    pshape = bayestopt_.pshape;
    p6 = bayestopt_.p6;
    p7 = bayestopt_.p7;
    p3 = bayestopt_.p3;
    p4 = bayestopt_.p4;
    npar = length(p6);
    pdraw = zeros(npar,1);
    lbcum = zeros(npar,1);
    ubcum = ones(npar,1);
    [junk1,junk2,junk3,lb,ub,junk4] = set_prior(estim_params_,M_,options_); %Prepare bounds
    if ~isempty(bayestopt_) && any(bayestopt_.pshape > 0)
        % Set prior bounds
        bounds = prior_bounds(bayestopt_, options_.prior_trunc);
        bounds.lb = max(bounds.lb,lb);
        bounds.ub = min(bounds.ub,ub);
    else  % estimated parameters but no declared priors
        % No priors are declared so Dynare will estimate the model by
        % maximum likelihood with inequality constraints for the parameters.
        bounds.lb = lb;
        bounds.ub = ub;
    end
    % set bounds for cumulative probabilities
    for i = 1:npar
      switch pshape(i)
        case 5% Uniform prior.
          p4(i) = min(p4(i),bounds.ub(i));
          p3(i) = max(p3(i),bounds.lb(i));
        case 3% Gaussian prior.
          lbcum(i) = 0.5 * erfc(-(bounds.lb(i)-p6(i))/p7(i) ./ sqrt(2));
          ubcum(i) = 0.5 * erfc(-(bounds.ub(i)-p6(i))/p7(i) ./ sqrt(2));
        case 2% Gamma prior.
          lbcum(i) = gamcdf(bounds.lb(i)-p3(i),p6(i),p7(i));
          ubcum(i) = gamcdf(bounds.ub(i)-p3(i),p6(i),p7(i));
        case 1% Beta distribution (TODO: generalized beta distribution)
          lbcum(i) = betainc((bounds.lb(i)-p3(i))./(p4(i)-p3(i)),p6(i),p7(i));
          ubcum(i) = betainc((bounds.ub(i)-p3(i))./(p4(i)-p3(i)),p6(i),p7(i));
        case 4% INV-GAMMA1 distribution
          % TO BE CHECKED
          lbcum(i) = gamcdf(1/(bounds.ub(i)-p3(i))^2,p7(i)/2,2/p6(i));
          ubcum(i) = gamcdf(1/(bounds.lb(i)-p3(i))^2,p7(i)/2,2/p6(i));
        case 6% INV-GAMMA2 distribution
          % TO BE CHECKED
          lbcum(i) = gamcdf(1/(bounds.ub(i)-p3(i)),p7(i)/2,2/p6(i));
          ubcum(i) = gamcdf(1/(bounds.lb(i)-p3(i)),p7(i)/2,2/p6(i));
        case 8
          lbcum(i) = weibcdf(bounds.lb(i)-p3(i),p6(i),p7(i));
          ubcum(i) = weibcdf(bounds.ub(i)-p3(i),p6(i),p7(i));
        otherwise
          % Nothing to do here.
      end
    end
    return
end


for i = 1:npar
    rdraw(:,i) = rdraw(:,i).*(ubcum(i)-lbcum(i))+lbcum(i);
    switch pshape(i)
      case 5% Uniform prior.
        pdraw(:,i) = rdraw(:,i)*(p4(i)-p3(i)) + p3(i);
      case 3% Gaussian prior.
        pdraw(:,i) = norminv(rdraw(:,i),p6(i),p7(i));
      case 2% Gamma prior.
        pdraw(:,i) = gaminv(rdraw(:,i),p6(i),p7(i))+p3(i);
      case 1% Beta distribution (TODO: generalized beta distribution)
        pdraw(:,i) = betainv(rdraw(:,i),p6(i),p7(i))*(p4(i)-p3(i))+p3(i);
      case 4% INV-GAMMA1 distribution 
        % TO BE CHECKED
        pdraw(:,i) =  sqrt(1./gaminv(rdraw(:,i),p7(i)/2,2/p6(i)))+p3(i);
      case 6% INV-GAMMA2 distribution  
        % TO BE CHECKED
        pdraw(:,i) =  1./gaminv(rdraw(:,i),p7(i)/2,2/p6(i))+p3(i);
      case 8
        pdraw(:,i) =  wblinv(rdraw(:,i),p6(i),p7(i))+p3(i);
      otherwise
        % Nothing to do here.
    end
end

  
