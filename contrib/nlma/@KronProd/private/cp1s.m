function out=cp1s(in)
% Outputs a matrix of ones the same size as IN


if islogical(in)

 out=true(size(in));

else,

 out=ones(size(in),class(in));

% else
% 
%  out=ones(size(in));

end
