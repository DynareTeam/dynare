function [xparams, logpost, options_]=metropolis_draw(init,options_,estim_params_,M_)
% function [xparams, logpost]=metropolis_draw(init) 
% Builds draws from metropolis
%
% INPUTS:
%   init:              scalar equal to 
%                      1: first call to store the required information 
%                           on files, lines, and chains to be read 
%                           in persistent variables to make them available 
%                           for future calls
%                      0: load a parameter draw
% Additional Inputs required for initialization
%   options_           [structure]     Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   M_                 [structure]     Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%   estim_params_      [structure]     Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%
% OUTPUTS:
%   xparams:           if init==0: vector of estimated parameters
%                      if init==1: error flaog
%   logpost:           log of posterior density
%   options_:          [structure]     Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   
% SPECIAL REQUIREMENTS
%
%   Requires CutSample to be run before in order to set up mh_history-file

% Copyright (C) 2003-2015 Dynare Team
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

persistent mh_nblck NumberOfDraws BaseName FirstLine FirstMhFile MAX_nruns

xparams = 0;
logpost = 0;

if init
    %get number of parameters
    nvx  = estim_params_.nvx;
    nvn  = estim_params_.nvn;
    ncx  = estim_params_.ncx;
    ncn  = estim_params_.ncn;
    np   = estim_params_.np ;
    npar = nvx+nvn+ncx+ncn+np;
    %get path of metropolis files
    MetropolisFolder = CheckPath('metropolis',M_.dname);
    FileName = M_.fname;
    BaseName = [MetropolisFolder filesep FileName];
    %load mh_history-file with info on what to load
    load_last_mh_history_file(MetropolisFolder, FileName);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine; 
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); 
    LastMhFile = TotalNumberOfMhFiles; 
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8); %number of parameters plus posterior plus ?
    mh_nblck = options_.mh_nblck;
    % set sub_draws option if empty
    if isempty(options_.sub_draws)
        options_.sub_draws = min(options_.posterior_max_subsample_draws, round(NumberOfDraws*mh_nblck));
    else
        if options_.sub_draws>NumberOfDraws*mh_nblck
            skipline()
            disp(['Estimation::mcmc: The value of option sub_draws (' num2str(options_.sub_draws) ') is greater than the number of available draws in the MCMC (' num2str(NumberOfDraws*mh_nblck) ')!'])
            disp('Estimation::mcmc: You can either change the value of sub_draws, reduce the value of mh_drop, or run another mcmc (with the load_mh_file option).')
            skipline()
            xparams = 1; % xparams is interpreted as an error flag
        end
    end
    return
else %not initialization, return one draw
    %get random draw from random chain
    ChainNumber = ceil(rand*mh_nblck);
    DrawNumber  = ceil(rand*NumberOfDraws);

    if DrawNumber <= MAX_nruns-FirstLine+1 %draw in first file, needs to account for first line
        MhFilNumber = FirstMhFile;
        MhLine = FirstLine+DrawNumber-1;
    else %draw in other file 
        DrawNumber  = DrawNumber-(MAX_nruns-FirstLine+1);
        MhFilNumber = FirstMhFile+ceil(DrawNumber/MAX_nruns); 
        MhLine = DrawNumber-(MhFilNumber-FirstMhFile-1)*MAX_nruns;
    end
    %load parameters and posterior
    load( [ BaseName '_mh' int2str(MhFilNumber) '_blck' int2str(ChainNumber) '.mat' ],'x2','logpo2');
    xparams = x2(MhLine,:);
    logpost= logpo2(MhLine);
end