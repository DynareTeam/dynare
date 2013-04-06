function set_dynare_seed(a,b)
% Set seeds depending on matlab (octave) version. This routine is called in dynare_config and can be called by the 
% user in the mod file.
%    
% Copyright (C) 2010-2012 Dynare Team
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
global options_

if ~nargin
    error('set_dynare_seed:: I need at least one input argument!')
end

matlab_random_streams = ~(exist('OCTAVE_VERSION') || matlab_ver_less_than('7.7') || options_.parallel_info.isHybridMatlabOctave);

if matlab_random_streams% Use new matlab interface.
    if nargin==1
        if ischar(a) && strcmpi(a,'default')
            options_.DynareRandomStreams.algo = 'mt19937ar';
            options_.DynareRandomStreams.seed = 0;
            s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
            if matlab_ver_less_than('7.12')
                reset(RandStream.setDefaultStream(s));
            else
                reset(RandStream.setGlobalStream(s));
            end
            return
        end
        if ischar(a) && strcmpi(a,'reset')
            s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
            if matlab_ver_less_than('7.12')
                reset(RandStream.setDefaultStream(s));
            else
                reset(RandStream.setGlobalStream(s));
            end
            return
        end
        if ischar(a)
            error('set_dynare_seed:: something is wrong in the calling sequence!')
        end
        if ~ischar(a)
            options_.DynareRandomStreams.algo = 'mt19937ar';
            options_.DynareRandomStreams.seed = a;
            s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
            if matlab_ver_less_than('7.12')
                reset(RandStream.setDefaultStream(s));
            else
                reset(RandStream.setGlobalStream(s));
            end
            return
        end
    elseif nargin==2
        if ~ischar(a) || ~( strcmpi(a,'mcg16807') || ...
                            strcmpi(a,'mlfg6331_64') || ...
                            strcmpi(a,'mrg32k3a') || ...
                            strcmpi(a,'mt19937ar') || ...
                            strcmpi(a,'shr3cong') || ...
                            strcmpi(a,'swb2712') )
            disp('set_dynare_seed:: First argument must be string designing the uniform random number algorithm!')
            RandStream.list
            disp(' ')
            disp('set_dynare_seed:: Change the first input accordingly...')
            disp(' ')
            error(' ')
        end
        if ~isint(b)
            error('set_dynare_seed:: The second input argument must be an integer!')
        end
        options_.DynareRandomStreams.algo = a;
        options_.DynareRandomStreams.seed = b;
        s = RandStream(options_.DynareRandomStreams.algo,'Seed',options_.DynareRandomStreams.seed);
        if matlab_ver_less_than('7.12')
            reset(RandStream.setDefaultStream(s));
        else
            reset(RandStream.setGlobalStream(s));
        end
    end
else% Use old matlab interface.
    if nargin==1
        if ischar(a) && strcmpi(a,'default')
            if exist('OCTAVE_VERSION') || matlab_ver_less_than('7.4')
                options_.DynareRandomStreams.algo = 'state';
            else
                % Twister was introduced in MATLAB 7.4
                options_.DynareRandomStreams.algo = 'twister';
            end
            options_.DynareRandomStreams.seed = 0;
            rand(options_.DynareRandomStreams.algo,options_.DynareRandomStreams.seed);
            randn('state',options_.DynareRandomStreams.seed);
            return
        end
        if ischar(a) && strcmpi(a,'reset')
            rand(options_.DynareRandomStreams.algo,options_.DynareRandomStreams.seed);
            randn('state',options_.DynareRandomStreams.seed);
            return
        end
        if ~ischar(a) && isint(a)
            options_.DynareRandomStreams.seed = a;
            rand(options_.DynareRandomStreams.algo,options_.DynareRandomStreams.seed);
            randn('state',options_.DynareRandomStreams.seed);
        else
            error('set_dynare_seed:: Something is wrong in the calling sequence!')
        end
    else
        error('set_dynare_seed:: Cannot use more than one input argument with your version of Matlab/Octave!')
    end
end