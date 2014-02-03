function terms = strsplit(s, delimiter)
%function terms = strsplit(s, delimiter)
%STRSPLIT Splits a string into multiple terms
%
%   terms = strsplit(s)
%       splits the string s into multiple terms that are separated by
%       white spaces (white spaces also include tab and newline).
%
%       The extracted terms are returned in form of a cell array of
%       strings.
%
%   terms = strsplit(s, delimiter)
%       splits the string s into multiple terms that are separated by
%       the specified delimiter. 
%   
%   Remarks
%   -------
%       - Note that the spaces surrounding the delimiter are considered
%         part of the delimiter, and thus removed from the extracted
%         terms.
%
%       - If there are two consecutive non-whitespace delimiters, it is
%         regarded that there is an empty-string term between them.         
%
%   Examples
%   --------
%       % extract the words delimited by white spaces
%       ts = strsplit('I am using MATLAB');
%       ts <- {'I', 'am', 'using', 'MATLAB'}
%
%       % split operands delimited by '+'
%       ts = strsplit('1+2+3+4', '+');
%       ts <- {'1', '2', '3', '4'}
%
%       % It still works if there are spaces surrounding the delimiter
%       ts = strsplit('1 + 2 + 3 + 4', '+');
%       ts <- {'1', '2', '3', '4'}
%
%       % Consecutive delimiters results in empty terms
%       ts = strsplit('C,Java, C++ ,, Python, MATLAB', ',');
%       ts <- {'C', 'Java', 'C++', '', 'Python', 'MATLAB'}
%
%       % When no delimiter is presented, the entire string is considered
%       % as a single term
%       ts = strsplit('YouAndMe');
%       ts <- {'YouAndMe'}
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 9, 2008
%
% Copyright (c) 2009, Dahua Lin
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

%% parse and verify input arguments

assert(ischar(s) && ndims(s) == 2 && size(s,1) <= 1, ...
    'strsplit:invalidarg', ...
    'The first input argument should be a char string.');

if nargin < 2
    by_space = true;
else
    d = delimiter;
    assert(ischar(d) && ndims(d) == 2 && size(d,1) == 1 && ~isempty(d), ...
        'strsplit:invalidarg', ...
        'The delimiter should be a non-empty char string.');
    
    d = strtrim(d);
    by_space = isempty(d);
end
    
%% main

s = strtrim(s);

if by_space
    w = isspace(s);            
    if any(w)
        % decide the positions of terms        
        dw = diff(w);
        sp = [1, find(dw == -1) + 1];     % start positions of terms
        ep = [find(dw == 1), length(s)];  % end positions of terms
        
        % extract the terms        
        nt = numel(sp);
        terms = cell(1, nt);
        for i = 1 : nt
            terms{i} = s(sp(i):ep(i));
        end                
    else
        terms = {s};
    end
    
else    
    p = strfind(s, d);
    if ~isempty(p)        
        % extract the terms        
        nt = numel(p) + 1;
        terms = cell(1, nt);
        sp = 1;
        dl = length(delimiter);
        for i = 1 : nt-1
            terms{i} = strtrim(s(sp:p(i)-1));
            sp = p(i) + dl;
        end         
        terms{nt} = strtrim(s(sp:end));
    else
        terms = {s};
    end        
end

