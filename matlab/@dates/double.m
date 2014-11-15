function [B, C]  = double(A) % --*-- Unitary tests --*--

% Returns a vector of doubles with the fractional part corresponding
% to the subperiod. Used for plots and to store dates in a matrix.
%
% INPUTS 
%  o A     dates object.
%
% OUTPUTS  
%  o B     A.ndat*1 vector of doubles.
%  o C     integer scalar, the frequency (1, 4, 12 or 52).
%
% REMARKS 
%  Obviously the frequency is lost during the conversion.
    
% Copyright (C) 2013 Dynare Team
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

B = A.time(:,1)+(A.time(:,2)-1)/A.freq;
if nargout>1
    C = A.freq;
end

%@test:1
%$ % Define a dates object
%$ qq = dates('Q');
%$ B = qq(1950,1):qq(1951,1);
%$
%$ % Call the tested routine.
%$ try
%$     C = double(B);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Define expected results.
%$ E = [ones(4,1)*1950; 1951];
%$ E = E + [(transpose(1:4)-1)/4; 0];
%$ if t(1)
%$     t(2) = dassert(C,E);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a dates object
%$ qq = dates('Q');
%$
%$ % Call the tested routine.
%$ try
%$     C = NaN(2,1);
%$     C(1) = double(qq(1950,1));
%$     C(2) = double(qq(1950,2));
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Define expected results.
%$ E = ones(2,1)*1950;
%$ E = E + [0; .25];
%$ if t(1)
%$     t(2) = dassert(C,E);
%$ end
%$ T = all(t);
%@eof:2