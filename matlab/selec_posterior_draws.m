function SampleAddress = selec_posterior_draws(SampleSize,drsize)
% Selects a sample of draws from the posterior distribution and if nargin>1
% saves the draws in _pdraws mat files (metropolis folder). If drsize>0
% the dr structure, associated to the parameters, is also saved in _pdraws.
% This routine is more efficient than metropolis_draw.m because here an
% _mh file cannot be opened twice.
%
% INPUTS
%   o SampleSize     [integer]  Size of the sample to build.
%   o drsize         [double]   structure dr is drsize megaoctets.
%
% OUTPUTS
%   o SampleAddress  [integer]  A (SampleSize*4) array, each line specifies the
%                               location of a posterior draw:
%                                  Column 2 --> Chain number
%                                  Column 3 --> (mh) File number
%                                  Column 4 --> (mh) line number
%
% SPECIAL REQUIREMENTS
%   None.
%

% Copyright (C) 2006-2011 Dynare Team
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

global M_ options_ estim_params_ oo_

% Number of parameters:
npar = estim_params_.nvx;
npar = npar + estim_params_.nvn;
npar = npar + estim_params_.ncx;
npar = npar + estim_params_.ncn;
npar = npar + estim_params_.np;

% Select one task:
switch nargin
  case 1
    info = 0;
  case 2
    MAX_mega_bytes = 10;% Should be an option...
    if drsize>0
        info=2;
    else
        info=1;
    end
    drawsize = drsize+npar*8/1048576;
  otherwise
    error(['selec_posterior_draws:: Unexpected number of input arguments!'])
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

% Get informations about the mcmc:
load_last_mh_history_file(MetropolisFolder, ModelName);
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
LastMhFile = TotalNumberOfMhFiles;
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
mh_nblck = options_.mh_nblck;

% Randomly select draws in the posterior distribution:
SampleAddress = zeros(SampleSize,4);
for i = 1:SampleSize
    ChainNumber = ceil(rand*mh_nblck);
    DrawNumber  = ceil(rand*NumberOfDraws);
    SampleAddress(i,1) = DrawNumber;
    SampleAddress(i,2) = ChainNumber;
    if DrawNumber <= MAX_nruns-FirstLine+1
        MhFileNumber = FirstMhFile;
        MhLineNumber = FirstLine+DrawNumber-1;
    else
        DrawNumber  = DrawNumber-(MAX_nruns-FirstLine+1);
        MhFileNumber = FirstMhFile+ceil(DrawNumber/MAX_nruns);
        MhLineNumber = DrawNumber-(MhFileNumber-FirstMhFile-1)*MAX_nruns;
    end
    SampleAddress(i,3) = MhFileNumber;
    SampleAddress(i,4) = MhLineNumber;
end
SampleAddress = sortrows(SampleAddress,[3 2]);

% Selected draws in the posterior distribution, and if drsize>0
% reduced form solutions, are saved on disk.
if info
    %delete old stale files before creating new ones
    delete_stale_file([BaseName '_posterior_draws*.mat'])
    if  SampleSize*drawsize <= MAX_mega_bytes% The posterior draws are saved in one file.
        pdraws = cell(SampleSize,info);
        old_mhfile = 0;
        old_mhblck = 0;
        for i = 1:SampleSize
            mhfile = SampleAddress(i,3);
            mhblck = SampleAddress(i,2);
            if (mhfile ~= old_mhfile) || (mhblck ~= old_mhblck)
                load([BaseName '_mh' num2str(mhfile) '_blck' num2str(mhblck) '.mat'],'x2')
            end
            pdraws(i,1) = {x2(SampleAddress(i,4),:)};
            if info-1
                set_parameters(pdraws{i,1});
                [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
                pdraws(i,2) = { dr };
            end
            old_mhfile = mhfile;
            old_mhblck = mhblck;
        end
        clear('x2')
        save([BaseName '_posterior_draws1.mat'],'pdraws')
    else% The posterior draws are saved in xx files.
        NumberOfDrawsPerFile = fix(MAX_mega_bytes/drawsize);
        NumberOfFiles = ceil(SampleSize*drawsize/MAX_mega_bytes);
        NumberOfLines = SampleSize - (NumberOfFiles-1)*NumberOfDrawsPerFile;
        linee = 0;
        fnum  = 1;
        pdraws = cell(NumberOfDrawsPerFile,info);
        old_mhfile = 0;
        old_mhblck = 0;
        for i=1:SampleSize
            linee = linee+1;
            mhfile = SampleAddress(i,3);
            mhblck = SampleAddress(i,2);
            if (mhfile ~= old_mhfile) || (mhblck ~= old_mhblck)
                load([BaseName '_mh' num2str(mhfile) '_blck' num2str(mhblck) '.mat'],'x2')
            end
            pdraws(linee,1) = {x2(SampleAddress(i,4),:)};
            if info-1
                set_parameters(pdraws{linee,1});
                [dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_);
                pdraws(linee,2) = { dr };
            end
            old_mhfile = mhfile;
            old_mhblck = mhblck;
            if fnum < NumberOfFiles && linee == NumberOfDrawsPerFile
                linee = 0;
                save([BaseName '_posterior_draws' num2str(fnum) '.mat'],'pdraws')
                fnum = fnum+1;
                if fnum < NumberOfFiles
                    pdraws = cell(NumberOfDrawsPerFile,info);
                else
                    pdraws = cell(NumberOfLines,info);
                end
            end
        end
        save([BaseName '_posterior_draws' num2str(fnum) '.mat'],'pdraws')
    end
end