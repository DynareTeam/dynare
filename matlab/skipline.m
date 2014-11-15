function skipline(n, fid)
% function skipline(n, fid)
% This function prints n newlines to fid
%
% INPUTS 
%
%   n        [integer]    Number of newlines to print
%   fid      [integer]    file id returned by fopen
%    
% OUTPUTS 
%   None
%

% Copyright (C) 2013-2014 Dynare Team
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

if nargin < 2
	fid = 1;
	if nargin < 1
		n = 1;
	end
end

for i=1:n
	fprintf(fid,'\n');
end
end