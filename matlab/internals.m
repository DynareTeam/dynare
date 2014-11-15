function internals(flag, varargin)

%@info:
%! @deftypefn {Function File} internals (@var{flag},@var{a},@var{b}, ...)
%! @anchor{internals}
%! @sp 1
%! This command provides internal documentation and unitary tests for the matlab routines.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item flag
%! Mandatory argument: --doc (for displaying internal documentation) or --test (for performing unitary tests).
%! @item b
%! Name of the routine to be tested or for which internal documentation is needed.
%! @item c
%! Name of the routine to be tested.
%! @item d
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! None.
%! @sp 2
%! @strong{Examples}
%! @sp 1
%! The following instruction:
%! @sp 1
%! @example
%! internals --info particle/local_state_iteration
%! @end example
%! will display the internal documentation of the routine local_state_iteration located in the particle subfolder of the matlab directory.
%! @sp 1
%! The following instruction:
%! @sp 1
%! @example
%! internals --test particle/local_state_iteration
%! @end example
%! will execute the unitary tests associated the routine local_state_iteration.
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! [1] It is not possible to display the internal documentation of more than one routine.
%! @sp 1
%! [2] It is possible to perform unitary tests on a list of routines.
%! @sp 1
%! [3] For displaying the internal documentation, matlab calls texinfo which has to be installed.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! None.
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{utilities/tests/dynTest} @ref{utilities/doc/dynInfo}
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
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

more off

if strcmpi(flag,'--test')
    if nargin>1
        dynare_path = dynare_config([],0);
        number_of_matlab_routines = length(varargin);
        for i=1:number_of_matlab_routines
            dtest(varargin{i},[dynare_path '..' filesep 'tests']);
        end
    else
        disp('You have to specify at least one Matlab routine after --test flag!')
    end
    return
end

if strcmpi(flag,'--load-mh-history') || strcmpi(flag,'--display-mh-history')
    switch length(varargin)
      case 3
        fname = varargin{1};
        if ~isequal(varargin{2},'in')
            error('internals:: Calling sequence must be of the form: internals --load-mh-history fname in dname')
        end
        dname = varargin{3};
      case 1
        fname = varargin{1};
        dname = varargin{1};
      otherwise
        error('internals:: Wrong calling sequence! You should read the manual...')
    end
    o = load_last_mh_history_file([dname filesep 'metropolis'],fname);
    if strcmpi(flag,'--load-mh-history')
        assignin('caller','mcmc_informations',o);
    else
        oo = load_first_mh_history_file([dname filesep 'metropolis'],fname);
        local = load([fname '_results'],'bayestopt_');
        names = local.bayestopt_.name; %evalin('base','bayestopt_.name');
        str = ['MCMC set-up for ' fname ' mod file'];
        ltr = length(str);
        skipline()
        disp(repmat('=',1,ltr))
        disp(str)
        disp(repmat('=',1,ltr))
        skipline(2)
        oar = compute_overall_acceptance_ratio([dname filesep 'metropolis'],fname);
        for b=1:o.Nblck
            str = ['MCMC chain number ' num2str(b) ':'];
            ltr = length(str);
            disp(str);
            disp(repmat('-',1,ltr));
            skipline()
            disp([' o Number of MCMC files is ' num2str(sum(o.MhDraws(:,2)))]);
            disp([' o Number of draws per chain is ' num2str(sum(o.MhDraws(:,1)))]);
            disp([' o Acceptance ratio in the current chain is ' num2str(oar(b)*100,'%5.2f') '%']);
            disp([' o Initial value of the posterior kernel is: ' num2str(oo.InitialLogPost(b),'%10.5f')])
            disp([' o Last value of the posterior kernel is: ' num2str(o.LastLogPost(b),'%10.5f')])
            disp([' o State of the chain:'])
            skipline()
            d1 = num2str(transpose(oo.InitialParameters(b,:)),'%10.5f\n');
            d2 = num2str(transpose(o.LastParameters(b,:)),'%10.5f\n');
            d1s = size(d1,2);
            d2s = size(d2,2);
            c0 = repmat('   ',length(names)+2,1);
            c1 = char(' ', repmat('+',1,size(char(names),2)), char(names));
            s1 = char(' || ','++++',repmat(' || ', length(names),1));
            t1 = repmat(' ',1,d1s);
            if d1s<=7
                t1 = 'Initial';
            else
                diff = d1s-7;
                if isequal(mod(diff,2),0)
                    start = diff/2+1;
                else
                    start = (diff-1)/2+1;
                end
                t1(start:start+6) = 'Initial';
            end
            c2 = char(t1,repmat('+',1,size(d1,2)),d1);
            s2 = char(' | ','+++',repmat(' | ', length(names),1));
            t2 = repmat(' ',1,d2s);
            if d2s<=7
                t2 = 'Current';
            else
                diff = d2s-7;
                if isequal(mod(diff,2),0)
                    start = diff/2+1;
                else
                    start = (diff-1)/2+1;
                end
                t2(start:start+6) = 'Current';
            end
            c3 = char(t2,repmat('+',1,size(d2,2)), d2);
            disp([c0, c1, s1, c2, s2, c3]);
            skipline()
        end
    end
    return
end

if strcmpi(flag,'--info')
    if nargin==2
        dynare_config([],0);
        dynInfo(varargin{1})
    else
        if nargin<2
            disp('You have to specify a Matlab routine after --info flag!')
        else
            disp('I can only show internal documentation for one Matlab routine!')
        end
    end
    return
end

disp('You should read the manual...')
