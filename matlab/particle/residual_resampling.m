function return_resample = residual_resampling(particles,weights,noise)
% Resamples particles.

%@info:
%! @deftypefn {Function File} {@var{indx} =} residual_resampling (@var{weights})
%! @anchor{particle/residual_resampling}
%! @sp 1
%! Resamples particles.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item weights
%! n*1 vector of doubles, particles' weights.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item indx
%! n*1 vector of intergers, indices.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{particle/resample}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2013 Dynare Team
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

% AUTHOR(S) frederic DOT karame AT univ DASH evry DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

% What is the number of particles?
number_of_particles = length(weights);

switch length(noise)
  case 1
    kitagawa_resampling = 1;
  case number_of_particles
    kitagawa_resampling = 0;
  otherwise
    error(['particle::resampling: Unknown method! The size of the second argument (' inputname(noise) ') is wrong.'])
end

% Set vectors of indices.
jndx = 1:number_of_particles;
indx = zeros(1,number_of_particles);

% Multiply the weights by the number of particles.
WEIGHTS = number_of_particles*weights;

% Compute the integer part of the normalized weights.
iWEIGHTS = fix(WEIGHTS);

% Compute the number of resample
number_of_trials = number_of_particles-sum(iWEIGHTS);

if number_of_trials
  WEIGHTS = (WEIGHTS-iWEIGHTS)/number_of_trials;
  EmpiricalCDF = cumsum(WEIGHTS);
  if kitagawa_resampling
    u = (transpose(1:number_of_trials)-1+noise(:))/number_of_trials;  
  else 
    u = fliplr(cumprod(noise(1:number_of_trials).^(1./(number_of_trials:-1:1))));
  end
  j=1;
  for i=1:number_of_trials
    while (u(i)>EmpiricalCDF(j))
      j=j+1;
    end
    iWEIGHTS(j)=iWEIGHTS(j)+1;
    if kitagawa_resampling==0
       j=1; 
    end
  end
end

k=1;
for i=1:number_of_particles
  if (iWEIGHTS(i)>0)
    for j=k:k+iWEIGHTS(i)-1
      indx(j) = jndx(i);
    end
  end
  k = k + iWEIGHTS(i);
end

if particles==0
    return_resample = indx ;
else
    return_resample = particles(indx,:) ;
end
%@test:1
%$ % Define the weights
%$ weights = randn(2000,1).^2;
%$ weights = weights/sum(weights);
%$ % Initialize t.
%$ t = ones(1,1);
%$
%$ try
%$     indx1 = residual_resampling(weights);
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define the weights
%$ weights = exp(randn(2000,1));
%$ weights = weights/sum(weights);
%$ % Initialize t.
%$ t = ones(1,1);
%$
%$ try
%$     indx1 = residual_resampling(weights);
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ T = all(t);
%@eof:2