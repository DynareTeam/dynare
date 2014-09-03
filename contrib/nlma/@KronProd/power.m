function L=power(L,R)
%POWER method for KronProd.
%
% out=power(L,R) or L.^R
%
%L must be a KronProd and R must be a scalar power.
%
%Algorithm: the exponentation is applied operand-by-operand to the elements
%of L.opset. Note that for non-integer R, and unless all the operand data 
%in L is non-negative, the result of full(L.^R) may be different
%from full(L).^R, due to the non-uniqueness of the latter. 

if ~isa(L,'KronProd') |  ~isnum(R)

 error('Bad input types');

end


L.scalarcoeff=full(double(L.scalarcoeff.^R));

for ii=1:length(L.opset)

 L.opset{ii}=L.opset{ii}.^R;

end
