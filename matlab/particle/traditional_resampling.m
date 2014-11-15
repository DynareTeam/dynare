function return_resample = traditional_resampling(particles,weights,noise)
% Resamples particles.

%@info:
%! @deftypefn {Function File} {@var{indx} =} traditional_resampling (@var{weights},@var{noise})
%! @anchor{particle/traditional_resampling}
%! @sp 1
%! Resamples particles (Resampling à la Kitagawa or stratified resampling).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item weights
%! n*1 vector of doubles, particles' weights.
%! @item noise
%! n*1 vector of doubles sampled from a [0,1] uniform distribution (stratified resampling) or scalar double
%! sampled from a [0,1] uniform distribution (Kitagawa resampling).
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

% Initialize the returned argument.
indx = ones(number_of_particles,1);

% Select method.
switch length(noise)
  case 1
    kitagawa_resampling = 1;
  case number_of_particles
    kitagawa_resampling = 0;
  otherwise
    error(['particle::resampling: Unknown method! The size of the second argument (' inputname(3) ') is wrong.'])
end

% Get the empirical  CDF.
c = cumsum(weights);

% Draw a starting point.
if kitagawa_resampling
    randvec = (transpose(1:number_of_particles)-1+noise(:))/number_of_particles ;
else 
    randvec = fliplr(cumprod(noise.^(1./(number_of_particles:-1:1))));
end

% Start at the bottom of the CDF
if kitagawa_resampling
    j = 1;
    for i=1:number_of_particles
        while (randvec(i)>c(j))
            j = j+1;
        end
        indx(i) = j;
    end
else
    for i=1:number_of_particles
        indx(i) = sum(randvec(i)>c);
    end
    % Matlab's indices start at 1...
    indx = indx+1;
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
%$ t = ones(2,1);
%$
%$ % First, try the stratified resampling.
%$ try
%$     indx1 = traditional_resampling(weights,rand(2000,1));
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Second, try the Kitagawa resampling.
%$ try
%$     indx2 = traditional_resampling(weights,rand);
%$ catch
%$     t(2) = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define the weights
%$ weights = exp(randn(2000,1));
%$ weights = weights/sum(weights);
%$ % Initialize t.
%$ t = ones(2,1);
%$
%$ % First, try the stratified resampling.
%$ try
%$     indx1 = traditional_resampling(weights,rand(2000,1));
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Second, try the Kitagawa resampling.
%$ try
%$     indx2 = traditional_resampling(weights,rand);
%$ catch
%$     t(2) = 0;
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ % Set the number of particles.
%$ number_of_particles = 20000;
%$
%$ show_plot = 0;
%$ show_time = 1;
%$
%$ % Define the weights
%$ weights = randn(number_of_particles,1).^2;
%$ weights = weights/sum(weights);
%$
%$ % Compute the empirical CDF
%$ c = cumsum(weights);
%$
%$ % Stratified resampling.
%$ noise  = rand(number_of_particles,1);
%$
%$ if show_time
%$     disp('Stratified resampling timing:')
%$     tic
%$ end
%$
%$ indx1  = traditional_resampling(weights,noise);
%$
%$ if show_time
%$     toc
%$     tic
%$ end
%$
%$ indx1_ = zeros(number_of_particles,1);
%$ randvec = (transpose(1:number_of_particles)-1+noise)/number_of_particles;
%$ for i=1:number_of_particles
%$     j = 1;
%$     while (randvec(i)>c(j))
%$         j = j + 1;
%$     end
%$     indx1_(i) = j;
%$ end
%$
%$ if show_time
%$     toc
%$ end
%$
%$ % Kitagawa's resampling.
%$ noise  = rand;
%$
%$ if show_time
%$     disp('Kitagawa''s resampling timing:')
%$     tic
%$ end
%$
%$ indx2  = traditional_resampling(weights,noise);
%$
%$ if show_time
%$     toc
%$     tic
%$ end
%$
%$ indx2_ = zeros(number_of_particles,1);
%$ randvec = (transpose(1:number_of_particles)-1+noise)/number_of_particles;
%$ j = 1;
%$ for i=1:number_of_particles
%$     while (randvec(i)>c(j))
%$         j = j + 1;
%$     end
%$     indx2_(i) = j;
%$ end
%$
%$ if show_time
%$     toc
%$ end
%$
%$ % REMARK
%$ % Note that the alternative code used in this test is sensibly faster than the code proposed
%$ % in the routine for the resampling scheme à la Kitagawa...
%$
%$ if show_plot
%$     plot(randvec,c,'-r'), hold on, plot([randvec(1),randvec(end)],[c(1),c(end)],'-k'), hold off, axis tight, box on
%$ end
%$
%$ % Check results.
%$ t(1) = dassert(indx1,indx1_);
%$ t(2) = dassert(indx2,indx2_);
%$ T = all(t);
%@eof:3