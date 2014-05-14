function C=encell(M,dim)
%Splits M into cell array along the dimension DIM. If DIM is not specified
%and M is a vector, the split is component-wise.
%
%Supports all array types.

if nargin<2
 if isvec(M)
  if iscol(M)
   dim=1;
  else
   dim=2;
  end
 else
  error('Not enough argument info.');
 end
end

ss=size(M);
ee=cp1(ss);
ss=mat2cell(ss,1,ee);
ss{dim}=row1s(ss{dim});

C=mat2cell(M,ss{:});
