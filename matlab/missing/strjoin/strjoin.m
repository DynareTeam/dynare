function return_val = strjoin (cellstring, delimiter) % --*-- Unitary tests --*--
%return_val = strjoin (cellstring, delimiter)

% Copyright (C) 2007 Muthiah Annamalai
% Copyright (C) 2013-2015 Ben Abbott
% Copyright (C) 2015 Dynare Team
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

% -*- texinfo -*-
% @deftypefn  {Function File} {@var{str} =} strjoin (@var{cstr})
% @deftypefnx {Function File} {@var{str} =} strjoin (@var{cstr}, @var{delimiter})
% Join the elements of the cell string array, @var{cstr}, into a single
% string.
%
% If no @var{delimiter} is specified, the elements of @var{cstr} are
% separated by a space.
%
% If @var{delimiter} is specified as a string, the cell string array is
% joined using the string.  Escape sequences are supported.
%
% If @var{delimiter} is a cell string array whose length is one less than
% @var{cstr}, then the elements of @var{cstr} are joined by interleaving the
% cell string elements of @var{delimiter}.  Escape sequences are not
% supported.
%
% @example
% @group
% strjoin (@{'Octave','Scilab','Lush','Yorick'@}, '*')
%       @result{} 'Octave*Scilab*Lush*Yorick'
% @end group
% @end example
% @seealso{strsplit}
% @end deftypefn

% Author: Muthiah Annamalai <muthiah.annamalai@uta.edu>
% Author: Ben Abbott <bpabbott@mac.com>


if (nargin == 1)
    delimiter = ' ';
elseif (nargin < 1 || nargin > 2)
    error('strjoin must be called with one or two arguments');
elseif ~(iscellstr (cellstring) && (ischar (delimiter) || iscellstr (delimiter)))
    error('The input arguments do not have the required types. See the help.');
end

if (numel(cellstring) == 1)
    return_val = cellstring{1};
    return;
end

if (ischar (delimiter))
    delimiter = sprintf(delimiter);
    delimiter = {delimiter};
end

num = numel(cellstring);
if (numel(delimiter) == 1 && num > 1)
    delimiter = repmat (delimiter, 1, num);
    delimiter(end) = {''};
elseif (num > 0 && numel (delimiter) ~= num - 1)
    error ('strjoin:cellstring_delimiter_mismatch. The number of delimiters does not match the number of strings')
else
    delimiter(end+1) = {''};
end

return_val=[];
if (num == 0)
    return_val = '';
else
    for ii=1:size(cellstring,2)
        return_val = [return_val,cellstring{ii},delimiter{ii}];
    end
end

%@test:1
%!t(1)=dassert (strjoin ({'hello'}, '-'), 'hello');
%!t(2)=dassert (strjoin ({'hello', 'world'}), 'hello world');
%!t(3)=dassert (strjoin ({'Octave', 'Scilab', 'Lush', 'Yorick'},'*'),'Octave*Scilab*Lush*Yorick');
%!t(4)=dassert (strjoin ({'space', 'comma', 'dash', 'semicolon', 'done'},{' ', ',', '-', ';'}), 'space comma,dash-semicolon;done');
%!t(5)=dassert (strjoin ({'Octave','Scilab'},'\n'),sprintf('Octave\nScilab'));
%!t(6)=dassert (strjoin ({'Octave','Scilab'},{'\n'}), 'Octave\nScilab');
%!t(7)=dassert (strjoin ({},'foo'), '');
%$ T = all(t);
%@eof:1