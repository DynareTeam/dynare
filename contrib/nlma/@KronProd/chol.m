function M=chol(M,varargin)
%CHOL method for KronProd class
%
%The Cholesky decomposition of a Kronecker product is distributive across the
%operands of the product, when the operands are square, positive definite.
%
%For a KronProd object, M, the syntax
%
%  K = chol(M,varargin{:}) 
%
%tries to exploit this by applying
%
%    [k{i}] = chol@double(M.opset{i},varargin{:}) 
%
%to each individual operand M.opset{i} of the KronProd object, M. 
%
%The output KronProd object K has K.opset=k. 
%
%NOTE 1: The method will fail if any of the operands of the input KronProd M
%are not square, positive definite, or if it has a non-positive
%scalar coefficient, M.scalarcoeff.
%
%NOTE 2: Currently, only 1 output argument is returned, unlike for chol@double.


chol_opset=cell(1,M.numops);

for ii=1:M.numops

 chol_opset{ii}= chol(M.opset{ii},varargin{:});

end


opinds=M.opinds;

scalarcoeff=sqrt(M.scalarcoeff);

if ~isreal(scalarcoeff) || scalarcoeff<=0
 error 'Scalar coefficient must be non-negative.'
end

domainsizes=M.domainsizes;

M=KronProd(chol_opset,opinds,domainsizes,scalarcoeff);
