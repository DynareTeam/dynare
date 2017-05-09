function data  = load_m_file_data_legacy(datafile, varobs)

% Copyright (C) 2014 Dynare Team
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

cXDHdrXnqo5KwwVpTRuc6OprAW = datafile(1:end-2);
[pathtocXDHdrXnqo5KwwVpTRuc6OprAW,cXDHdrXnqo5KwwVpTRuc6OprAW,junk] = fileparts(cXDHdrXnqo5KwwVpTRuc6OprAW);

if ~isempty(pathtocXDHdrXnqo5KwwVpTRuc6OprAW)
    OvMuQsJgjwzYG5Pni0TzU8Acb2YBJva = pwd();
    cd(pathtocXDHdrXnqo5KwwVpTRuc6OprAW);
end

eval(cXDHdrXnqo5KwwVpTRuc6OprAW);

if ~isempty(pathtocXDHdrXnqo5KwwVpTRuc6OprAW)
    cd(OvMuQsJgjwzYG5Pni0TzU8Acb2YBJva);
end

try
    if size(eval(varobs{1,1}),2)>size(eval(varobs{1,1}),1)
        data = dseries(eval(cellofstring4eval(varobs,'transpose')),[],varobs);
    else
        data = dseries(eval(cellofstring4eval(varobs)),[],varobs);
    end
catch
    errmsg = sprintf('makedataset: Check that all the variables listed in varobs exist in %s, have the same number of observations, and are all column or row vectors.',datafile);
    error(errmsg)
end

function str = cellofstring4eval(A,transpose_option)
    n = length(A);
    str = '[';
    for i=1:n-1
        if nargin==2 && strcmp(transpose_option,'transpose')
            str = [str, A{i}, ''','];
        else
            str = [str, A{i}, ','];
        end
    end
    if nargin==2 && strcmp(transpose_option,'transpose')
        str = [str, A{n}, ''']'];
    else
        str = [str, A{n}, ']'];
    end