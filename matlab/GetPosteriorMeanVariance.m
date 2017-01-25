function [mean,variance] = GetPosteriorMeanVariance(M,drop)

% Copyright (C) 2012-2016 Dynare Team
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
    
    MetropolisFolder = CheckPath('metropolis',M.dname);
    FileName = M.fname;
    BaseName = [MetropolisFolder filesep FileName];
    load_last_mh_history_file(MetropolisFolder, FileName);
    NbrDraws = sum(record.MhDraws(:,1));
    NbrFiles = sum(record.MhDraws(:,2));
    NbrBlocks = record.Nblck;
    mean = 0;
    variance = 0;
    
    NbrKeptDraws = 0;
    for i=1:NbrBlocks
        NbrDrawsCurrentBlock = 0;
        for j=1:NbrFiles
            o = load([BaseName '_mh' int2str(j) '_blck' int2str(i),'.mat']);
            NbrDrawsCurrentFile = size(o.x2,1);
            if NbrDrawsCurrentBlock + NbrDrawsCurrentFile <= drop*NbrDraws
                NbrDrawsCurrentBlock = NbrDrawsCurrentBlock + NbrDrawsCurrentFile;
                continue
            elseif NbrDrawsCurrentBlock < drop*NbrDraws
                FirstDraw = ceil(drop*NbrDraws - NbrDrawsCurrentBlock + 1);
                x2 = o.x2(FirstDraw:end,:);
            else
                x2 = o.x2;
            end
            NbrKeptDrawsCurrentFile = size(x2,1);
            %recursively compute mean and variance
            mean = (NbrKeptDraws*mean + sum(x2)')/(NbrKeptDraws+NbrKeptDrawsCurrentFile);
            x2Demeaned = bsxfun(@minus,x2,mean');
            variance = (NbrKeptDraws*variance + x2Demeaned'*x2Demeaned)/(NbrKeptDraws+NbrKeptDrawsCurrentFile);
            NbrDrawsCurrentBlock = NbrDrawsCurrentBlock + NbrDrawsCurrentFile;
            NbrKeptDraws = NbrKeptDraws + NbrKeptDrawsCurrentFile;
        end
    end
