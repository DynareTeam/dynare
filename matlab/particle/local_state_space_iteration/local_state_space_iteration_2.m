function [y,y_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,a,b,c)%yhat_,ss)

%@info:
%! @deftypefn {Function File} {@var{y}, @var{y_} =} local_state_equation_2 (@var{yhat},@var{epsilon}, @var{ghx}, @var{ghu}, @var{constant}, @var{ghxx}, @var{ghuu}, @var{ghxu}, @var{yhat_}, @var{ss})
%! @anchor{particle/local_state_space_iteration_2}
%! @sp 1
%! Given the states (y) and structural innovations (epsilon), this routine computes the level of selected endogenous variables when the
%! model is approximated by an order two taylor expansion around the deterministic steady state. Depending on the number of input/output
%! argument the pruning algorithm advocated by C. Sims is used or not (this version should not be used if the selected endogenous variables
%! are not the states of the model).
%!
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item yhat
%! n*1 vector of doubles, initial condition, where n is the number of state variables.
%! @item epsilon
%! q*1 vector of doubles, structural innovations.
%! @item ghx
%! m*n matrix of doubles, restricted dr.ghx where we only consider the lines corresponding to a subset of endogenous variables.
%! @item ghu
%! m*q matrix of doubles, restricted dr.ghu where we only consider the lines corresponding to a subset of endogenous variables.
%! @item constant
%! m*1 vector of doubles, deterministic steady state plus second order correction for a subset of endogenous variables.
%! @item ghxx
%! m*n² matrix of doubles, restricted dr.ghxx where we only consider the lines corresponding to a subset of endogenous variables.
%! @item ghuu
%! m*q² matrix of doubles, restricted dr.ghuu where we only consider the lines corresponding to a subset of endogenous variables.
%! @item ghxu
%! m*(nq) matrix of doubles, subset of dr.ghxu where we only consider the lines corresponding to a subset of endogenous variables.
%! @item yhat_
%! n*1 vector of doubles, spurious states for the pruning version.
%! @item ss
%! n*1 vector of doubles, steady state for the states.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item y
%! m*1 vector of doubles, selected endogenous variables.
%! @item y_
%! m*1 vector of doubles, update of the latent variables needed for the pruning version (first order update).
%! @end table
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! [1] If the function has 10 input arguments then it must have 2 output arguments (pruning version).
%! @sp 1
%! [2] If the function has 08 input arguments then it must have 1 output argument.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2012 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr
%           frederic DOT karame AT univ DASH evry DOT fr

if nargin==9
    pruning = 0; numthreads = a;
    if nargout>1
        error('local_state_space_iteration_2:: Numbers of input and output argument are inconsistent!')
    end
elseif nargin==11
    pruning = 1; numthreads = c; yhat_ = a; ss = b;
    if nargout~=2
        error('local_state_space_iteration_2:: Numbers of input and output argument are inconsistent!')
    end
else
    error('local_state_space_iteration_2:: Wrong number of input arguments!')
end

number_of_threads = numthreads;

switch pruning
  case 0
    for i =1:size(yhat,2)
        y(:,i) = constant + ghx*yhat(:,i) + ghu*epsilon(:,i) ...
                 + A_times_B_kronecker_C(.5*ghxx,yhat(:,i),number_of_threads)  ...
                 + A_times_B_kronecker_C(.5*ghuu,epsilon(:,i),number_of_threads) ...
                 + A_times_B_kronecker_C(ghxu,yhat(:,i),epsilon(:,i),number_of_threads);
    end
  case 1
    for i =1:size(yhat,2)
        y(:,i) = constant + ghx*yhat(:,i) + ghu*epsilon(:,i) ...
                 + A_times_B_kronecker_C(.5*ghxx,yhat_(:,i),number_of_threads)  ...
                 + A_times_B_kronecker_C(.5*ghuu,epsilon(:,i),number_of_threads) ...
                 + A_times_B_kronecker_C(ghxu,yhat_(:,i),epsilon(:,i),number_of_threads);
    end
    y_ = ghx*yhat_+ghu*epsilon;
    y_ = bsxfun(@plus,y_,ss);
end

%@test:1
%$ n = 2;
%$ q = 3;
%$
%$ yhat = zeros(n,1);
%$ epsilon = zeros(q,1);
%$ ghx = rand(n,n);
%$ ghu = rand(n,q);
%$ constant = ones(n,1);
%$ ghxx = rand(n,n*n);
%$ ghuu = rand(n,q*q);
%$ ghxu = rand(n,n*q);
%$ yhat_ = zeros(n,1);
%$ ss = ones(n,1);
%$
%$ % Call the tested routine.
%$ for i=1:10
%$     y1 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1);
%$     [y2,y2_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,ss,1);
%$ end
%$
%$ % Check the results.
%$ t(1) = dassert(y1,ones(n,1));
%$ t(2) = dassert(y2,ones(n,1));
%$ t(3) = dassert(y2_,ones(n,1));
%$ T = all(t);
%@eof:1

%@test:2
%$ old_path = pwd;
%$ cd([fileparts(which('dynare')) '/../tests/']);
%$ dynare('dsge_base2');
%$ load dsge_base2;
%$ cd(old_path);
%$ dr = oo_.dr;
%$ clear('oo_','options_','M_');
%$ delete([fileparts(which('dynare')) '/../tests/dsge_base2.mat']);
%$ istates = dr.nstatic+(1:dr.npred);
%$ n = dr.npred;
%$ q = size(dr.ghu,2);
%$ yhat = zeros(n,1);
%$ epsilon = zeros(q,1);
%$ ghx = dr.ghx(istates,:);
%$ ghu = dr.ghu(istates,:);
%$ constant = dr.ys(istates,:)+dr.ghs2(istates,:);
%$ ghxx = dr.ghxx(istates,:);
%$ ghuu = dr.ghuu(istates,:);
%$ ghxu = dr.ghxu(istates,:);
%$ yhat_ = zeros(n,1);
%$ ss = dr.ys(istates,:);
%$
%$ t = ones(2,1);
%$
%$ % Call the tested routine.
%$ try
%$     y1 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1);
%$ catch
%$     t(1) = 0;
%$ end
%$ try
%$     [y2,y2_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,ss,1);
%$ catch
%$     t(2) = 0;
%$ end
%$
%$ % Check the results.
%$ T = all(t);
%@eof:2

%@test:3
%$ Bohrbug = 1; % A bug that manifests reliably under a possibly unknown but well-defined set of conditions.
%$ if ~Bohrbug
%$     n = 2;
%$     q = 3;
%$
%$     yhat = .01*randn(n,1);
%$     epsilon = .001*randn(q,1);
%$     ghx = rand(n,n);
%$     ghu = rand(n,q);
%$     constant = ones(n,1);
%$     ghxx = rand(n,n*n);
%$     ghuu = rand(n,q*q);
%$     ghxu = rand(n,n*q);
%$     yhat_ = zeros(n,1);
%$     ss = ones(n,1);
%$
%$     % Call the tested routine (mex version).
%$     y1a = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1);
%$     [y2a,y2a_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,ss,1);
%$
%$     % Call the tested routine (matlab version)
%$     path_to_mex = fileparts(which(['qmc_sequence.' mexext]));
%$     where_am_i_coming_from = pwd;
%$     cd(path_to_mex);
%$     tar('local_state_space_iteration_2.tar',['local_state_space_iteration_2.' mexext]);
%$     cd(where_am_i_coming_from);
%$     dynare_config([],0);
%$     y1b = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1);
%$     [y2b,y2b_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,ss,1);
%$     cd(path_to_mex);
%$     untar('local_state_space_iteration_2.tar');
%$     delete('local_state_space_iteration_2.tar');
%$     cd(where_am_i_coming_from);
%$     dynare_config([],0);
%$     % Check the results.
%$     t(1) = dassert(y1a,y1b);
%$     t(2) = dassert(y2a,y2b);
%$     t(3) = dassert(y2a_,y2b_);
%$     T = all(t);
%$ else
%$     t(1) = 1;
%$     T = all(t);
%$ end
%@eof:3


%@test:4
%$ % TIMING TEST (parallelization with openmp)
%$ old_path = pwd;
%$ cd([fileparts(which('dynare')) '/../tests/']);
%$ dynare('dsge_base2');
%$ load dsge_base2;
%$ cd(old_path);
%$ dr = oo_.dr;
%$ clear('oo_','options_','M_');
%$ delete([fileparts(which('dynare')) '/../tests/dsge_base2.mat']);
%$ istates = dr.nstatic+(1:dr.npred);
%$ n = dr.npred;
%$ q = size(dr.ghu,2);
%$ yhat = zeros(n,10000000);
%$ epsilon = zeros(q,10000000);
%$ ghx = dr.ghx(istates,:);
%$ ghu = dr.ghu(istates,:);
%$ constant = dr.ys(istates,:)+dr.ghs2(istates,:);
%$ ghxx = dr.ghxx(istates,:);
%$ ghuu = dr.ghuu(istates,:);
%$ ghxu = dr.ghxu(istates,:);
%$ yhat_ = zeros(n,10000000);
%$ ss = dr.ys(istates,:);
%$
%$ t = NaN(4,1);
%$ tic, for i=1:10, y1 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1); end
%$ t1 = toc;
%$ tic, for i=1:10, y2 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,2); end
%$ t2 = toc;
%$ t(1) = dassert(y1,y2,1e-15); clear('y1');
%$ tic, for i=1:10, y3 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,3); end
%$ t3 = toc;
%$ t(2) = dassert(y2,y3,1e-15); clear('y2');
%$ tic, for i=1:10, y4 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,4); end
%$ t4 = toc;
%$ t(3) = dassert(y4,y3,1e-15); clear('y3','y4');
%$ t(4) = (t1>t2) && (t2>t3) && (t3>t4);
%$ if ~t(4)
%$     disp('Timmings:')
%$     [t1, t2, t3, t4]
%$ end
%$ % Check the results.
%$ T = all(t);
%@eof:4
