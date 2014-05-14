function y=enrow(x)
% Makes x into a row vector if it isn't already


y=x(:).';

%%%%%%%%%%%%%%%%%%%%
if ~isvec(x), error('Vector input required.'); end

if ~isrow(x), 
 y=x.'; 
else
y=x;
end
