function L=mpower(L,R)
%MPOWER method for KronProd.
%
% out=mpower(L,R) or out=L^R
%
%L must be a KronProd and R must be a scalar power
%
%Algorithm: the exponentation is applied operand-by-operand to the elements
%of L.opset.

if ~isa(L,'KronProd') |  ~isnum(R)

 error('Bad input types');

end

L.scalarcoeff=full(double(L.scalarcoeff^R));

for ii=1:length(L.opset)

 L.opset{ii}=L.opset{ii}^R;

end
