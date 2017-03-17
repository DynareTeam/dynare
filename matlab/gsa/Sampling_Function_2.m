function [Outmatrix, OutFact] = Sampling_Function_2(p, k, r, UB, LB, GroupMat)
%[Outmatrix, OutFact] = Sampling_Function_2(p, k, r, UB, LB, GroupMat)
%	Inputs: k (1,1)                      := number of factors examined or number of groups examined.
%                                           In case the groups are chosen the number of factors is stores in NumFact and
%                                           sizea becomes the number of created groups. 
%           NumFact (1,1)                := number of factors examined in the case when groups are chosen
%	    	r (1,1)                      := sample size  
%           p (1,1)                      := number of intervals considered in [0, 1]
%           UB(sizea,1)                  := Upper Bound for each factor 
%           LB(sizea,1)                  := Lower Bound for each factor 
%           GroupNumber(1,1)             := Number of groups (eventually 0)
%           GroupMat(NumFact,GroupNumber):= Matrix which describes the chosen groups. Each column represents a group and its elements 
%                                           are set to 1 in correspondence of the factors that belong to the fixed group. All
%                                           the other elements are zero.
%   Local Variables:  
%	    	sizeb (1,1)         := sizea+1
%           sizec (1,1)         := 1
%           randmult (sizea,1)  := vector of random +1 and -1  
%           perm_e(1,sizea)     := vector of sizea random permutated indeces    
%           fact(sizea)         := vector containing the factor varied within each traj
% 	        DDo(sizea,sizea)    := D*       in Morris, 1991   
%	        A(sizeb,sizea)      := Jk+1,k   in Morris, 1991
%	        B(sizeb,sizea)      := B        in Morris, 1991
%	        Po(sizea,sizea)     := P*       in Morris, 1991
%           Bo(sizeb,sizea)     := B*       in Morris, 1991
%	        Ao(sizeb,sizec)     := Jk+1,1   in Morris, 1991
%	        xo(sizec,sizea)     := x*       in Morris, 1991 (starting point for the trajectory)
%           In(sizeb,sizea)     := for each loop orientation matrix. It corresponds to a trajectory
%                                  of k step in the parameter space and it provides a single elementary
%                                  effect per factor 
%           MyInt()
%           Fact(sizea,1)       := for each loop vector indicating which factor or group of factors has been changed 
%                                  in each step of the trajectory
%           AuxMat(sizeb,sizea) := Delta*0.5*((2*B - A) * DD0 + A) in Morris, 1991. The AuxMat is used as in Morris design
%                                  for single factor analysis, while it constitutes an intermediate step for the group analysis.
%
%	Output: Outmatrix(sizeb*r, sizea) := for the entire sample size computed In(i,j) matrices
%           OutFact(sizea*r,1)        := for the entire sample size computed Fact(i,1) vectors
%           
%   Note: B0 is constructed as in Morris design when groups are not considered. When groups are considered the routine
%         follows the following steps:
%           1- Creation of P0 and DD0 matrices defined in Morris for the groups. This means that the dimensions of these
%              2 matrices are (GroupNumber,GroupNumber).
%           2- Creation of AuxMat matrix with (GroupNumber+1,GroupNumber) elements.
%           3- Definition of GroupB0 starting from AuxMat, GroupMat and P0.
%           4- The final B0 for groups is obtained as [ones(sizeb,1)*x0' + GroupB0]. The P0 permutation is present in GroupB0
%              and it's not necessary to permute the matrix (ones(sizeb,1)*x0') because it's already randomly created. 
%   Reference:
%   A. Saltelli, K. Chan, E.M. Scott, "Sensitivity Analysis" on page 68 ss
%
%   F. Campolongo, J. Cariboni, JRC - IPSC Ispra, Varese, IT
%   Last Update: 15 November 2005 by J.Cariboni
%
% Written by Jessica Cariboni and Francesca Campolongo
% Joint Research Centre, The European Commission,
%

% Copyright (C) 2005 European Commission
% Copyright (C) 2012 Dynare Team
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

% Parameters and initialisation of the output matrix
sizea = k;
Delta = p/(2*p-2);
%Delta = 1/3
NumFact = sizea;
GroupNumber = size(GroupMat,2);

if GroupNumber ~ 0;
    sizea = size(GroupMat,2);
end    

sizeb = sizea + 1;
sizec = 1;
Outmatrix = [];
OutFact = [];

% For each i generate a trajectory  
for i=1:r
    
    % Construct DD0 - OLD VERSION - it does not need communication toolbox
    % RAND(N,M) is an NXM matrix with random entries, chosen from a uniform distribution on the interval (0.0,1.0).
    % Note that DD0 tells if the factor have to be increased or ddecreased
    % by Delta.
    randmult = ones(k,1);           
    v = rand(k,1);                  
    randmult (find(v < 0.5))=-1;
    randmult = repmat(randmult,1,k);
    DD0 = randmult .* eye(k);
    
    % Construct DD0 - NEW VERSION - it needs communication toolbox
    % randsrc(m) generates an m-by-m matrix, each of whose entries independently takes the value -1 with probability 1/2,
    % and 1 with probability 1/2.
    % DD0 = randsrc(NumFact) .* eye(NumFact);      
    
    % Construct B (lower triangular)
    B = ones(sizeb,sizea);
    for j = 1:sizea
       B(1:j,j)=0;    
    end
    
    % Construct A0, A
    A0 = ones(sizeb,1);
    A = ones(sizeb,NumFact);

    % Construct the permutation matrix P0. In each column of P0 one randomly chosen element equals 1
    % while all the others equal zero. 
    % P0 tells the order in which order factors are changed in each
    % trajectory. P0 is created as it follows:
    % 1) All the elements of P0 are set equal to zero ==> P0 = zeros (sizea, sizea);
    % 2) The function randperm create a random permutation of integer 1:sizea, without repetitions ==> perm_e; 
    % 3) In each column of P0 the element indicated in perm_e is set equal to one.    
    % Note that P0 is then used reading it by rows. 
    P0 = zeros (sizea, sizea);
    perm_e = randperm(sizea);               % RANDPERM(n) is a random permutation of the integers from 1 to n.
    for j = 1:sizea
        P0(perm_e(j),j) = 1;    
    end    
    
    % When groups are present the random permutation is done only on B. The effect is the same since 
    % the added part (A0*x0') is completely random. 
    if GroupNumber ~ 0
        B = B * (GroupMat*P0')';
    end
    
    % Compute AuxMat both for single factors and groups analysis. For Single factors analysis
    % AuxMat is added to (A0*X0) and then permutated through P0. When groups are active AuxMat is
    % used to build GroupB0. AuxMat is created considering DD0. If the element on DD0 diagonal
    % is 1 then AuxMat will start with zero and add Delta. If the element on DD0 diagonal is -1 
    % then DD0 will start Delta and goes to zero.
    AuxMat = Delta*0.5*((2*B - A) * DD0 + A);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a --> Define the random vector x0 for the factors. Note that x0 takes value in the hypercube
    % [0,...,1-Delta]*[0,...,1-Delta]*[0,...,1-Delta]*[0,...,1-Delta] 
    MyInt = repmat([0:(1/(p-1)):(1-Delta)],NumFact,1);     % Construct all possible values of the factors
    
    % OLD VERSION - it needs communication toolbox
    % w = randint(NumFact,1,[1,size(MyInt,2)]);              
    
    % NEW VERSION - construct a version of random integers
    % 1) create a vector of random integers
    % 2) divide [0,1] into the needed steps
    % 3) check in which interval the random numbers fall
    % 4) generate the corresponding integer
    v = repmat(rand(NumFact,1),1,size(MyInt,2)+1);     % 1)
    IntUsed = repmat([0:1/size(MyInt,2):1],NumFact,1); % 2)
    DiffAuxVec = IntUsed - v;                          % 3)
    
    for ii = 1:size(DiffAuxVec,1)
        w(1,ii) = max(find(DiffAuxVec(ii,:)<0));       % 4)
    end
    x0 = MyInt(1,w)';                                  % Define x0    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % b --> Compute the matrix B*, here indicated as B0. Each row in B0 is a
    % trajectory for Morris Calculations. The dimension of B0 is (Numfactors+1,Numfactors) 
    if GroupNumber ~ 0
        B0 = (A0*x0' + AuxMat);
    else
        B0 = (A0*x0' + AuxMat)*P0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % c --> Compute values in the original intervals
    % B0 has values x(i,j) in [0, 1/(p -1), 2/(p -1), ... , 1].
    % To obtain values in the original intervals [LB, UB] we compute
    % LB(j) + x(i,j)*(UB(j)-LB(j))
    In = repmat(LB,1,sizeb)' + B0 .* repmat((UB-LB),1,sizeb)';

    % Create the Factor vector. Each component of this vector indicate which factor or group of factor
    % has been changed in each step of the trajectory.
    for j=1:sizea
        Fact(1,j) = find(P0(j,:));
    end
    Fact(1,sizea+1) = 0;
    
    Outmatrix = [Outmatrix; In];
    OutFact = [OutFact; Fact'];
    
end