function [SAmeas, OutMatrix] = Morris_Measure_Groups(NumFact, Sample, Output, p, Group)
% [SAmeas, OutMatrix] = Morris_Measure_Groups(NumFact, Sample, Output, p, Group)
%
% Given the Morris sample matrix, the output values and the group matrix compute the Morris measures
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% Group [NumFactor, NumGroups] := Matrix describing the groups. 
% Each column represents one group. 
% The element of each column are zero if the factor is not in the
% group. Otherwise it is 1.
%
% Sample := Matrix of the Morris sampled trajectories 
%
% Output := Matrix of the output(s) values in correspondence of each point
% of each trajectory
%
% k = Number of factors
% -------------------------------------------------------------------------
% OUTPUTS 
% OutMatrix (NumFactor*NumOutputs, 3)= [Mu*, Mu, StDev]
% for each output it gives the three measures of each factor
% -------------------------------------------------------------------------
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

if nargin==0,
  skipline()
  disp('[SAmeas, OutMatrix] = Morris_Measure_Groups(NumFact, Sample, Output, p, Group);')
  return
end

OutMatrix=[];
if nargin < 5, Group=[]; end

NumGroups = size(Group,2);
if nargin < 4 | isempty(p),
    p = 4;
end
Delt = p/(2*p-2);

if NumGroups ~ 0
    sizea = NumGroups;      % Number of groups
    GroupMat=Group;
    GroupMat = GroupMat';
else 
    sizea = NumFact; 
end
r=size(Sample,1)/(sizea+1);     % Number of trajectories

% For Each Output
for k=1:size(Output,2)
    
    OutValues=Output(:,k);
  
    % For each r trajectory
    for i=1:r
    
        % For each step j in the trajectory
        % Read the orientation matrix fact for the r-th sampling
        % Read the corresponding output values
        Single_Sample = Sample(i+(i-1)*sizea:i+(i-1)*sizea+sizea,:); 
        Single_OutValues = OutValues(i+(i-1)*sizea:i+(i-1)*sizea+sizea,:); 
        A = (Single_Sample(2:sizea+1,:)-Single_Sample(1:sizea,:))';
        Delta = A(find(A));

        % For each point of the fixed trajectory compute the values of the Morris function. The function
        % is partitioned in four parts, from order zero to order 4th.
        for j=1:sizea   % For each point in the trajectory i.e for each factor   
            % matrix of factor which changes
            if NumGroups ~ 0;
                AuxFind (:,1) = A(:,j);
%                 AuxFind(find(A(:,j)),1)=1;
%                 Pippo = sum((Group - repmat(AuxFind,1,NumGroups)),1);
%                 Change_factor(j,i) = find(Pippo==0);   
                Change_factor = find(abs(AuxFind)>1e-010); 
                % If we deal with groups we can only estimate the new mu*
                % measure since factors in the same groups can move in
                % opposite direction and the definition of the standard
                % Morris mu cannopt be applied. 
                % In the new version the elementary effect is defined with
                % the absolute value.
                %SAmeas(find(GroupMat(Change_factor(j,i),:)),i) = abs((Single_OutValues(j) - Single_OutValues(j+1) )/Delt); %(2/3));   
                SAmeas(i,Change_factor') = abs((Single_OutValues(j) - Single_OutValues(j+1) )/Delt);   
            else
                Change_factor(j,i) = find(Single_Sample(j+1,:)-Single_Sample(j,:));
                % If no groups --> we compute both the original and
                % modified measure
                if Delta(j) > 0                              %=> +Delta
                    SAmeas(Change_factor(j,i),i) = (Single_OutValues(j+1) - Single_OutValues(j) )/Delt; %(2/3);
                else                                         %=> -Delta
                    SAmeas(Change_factor(j,i),i) = (Single_OutValues(j) - Single_OutValues(j+1) )/Delt; %(2/3);
                end 
            end
        end   %for j=1:sizea
    
    end     %for i=1:r
   
    if NumGroups ~ 0
        SAmeas = SAmeas';
    end

    % Compute Mu AbsMu and StDev
    if any(any(isnan(SAmeas)))
      for j=1:NumFact,
        SAm = SAmeas(j,:);
        SAm = SAm(find(~isnan(SAm)));
        rr=length(SAm);
        AbsMu(j,1) = sum(abs(SAm),2)/rr;
      if NumGroups == 0
        Mu(j,1) = sum(SAm,2)/rr;
        StDev(j,1) = sum((SAm - repmat(Mu(j),1,rr)).^2/(rr*(rr-1)),2).^0.5;
      end
      end
    else
      AbsMu = sum(abs(SAmeas),2)/r;
      if NumGroups == 0
        Mu = sum(SAmeas,2)/r;
        StDev = sum((SAmeas - repmat(Mu,1,r)).^2/(r*(r-1)),2).^0.5;
      end
    end

    % Define the output Matrix - if we have groups we cannot define the old
    % measure mu, only mu* makes sense
    if NumGroups > 0
        OutMatrix = [OutMatrix; AbsMu];   
    else
        OutMatrix = [OutMatrix; AbsMu, Mu, StDev];   
    end
end     % For Each Output
