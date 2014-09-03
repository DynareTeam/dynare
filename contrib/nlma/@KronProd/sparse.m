function T=sparse(M)
%SPARSE - sparse() method for KronProd class. Convert to sparse numeric (double) form. 
%
%  T=sparse(M)
%
%Currently requires that M is a KronProd consisting of 
%numeric or logical matrix class operators only.
%
%The constituent operands are combined into one large matrix using kron().
%Scalar constituent operands are converted to a multiple of the sparse identity matrix
%when doing so.
%
%Note: this method is obviously not recommended when M.opset{} are large and dense matrices.




allnumeric=1;

for ii=1:length(M.opset)

 allnumeric=allnumeric&( isnumeric(M.opset{ii})  | islogical(M.opset{ii}) );

end




if ~allnumeric
 error('NUMERIC() method only defined when constituent operators are all numeric matrices');
end




T=M.scalarcoeff;
nn=M.maxdim;

while nn

  op=sparse(double(M.opset{M.opinds(nn)}));

  if isnum(op), op=speye(M.domainsizes(nn))*op; end


  T=kron(T,op);

  nn=nn-1;

end


