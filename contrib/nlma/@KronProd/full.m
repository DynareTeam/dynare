function M=full(M)
%FULL - full() method for KronProd class. Convert to full numeric (double) form. 
%
%  T=full(M)
%
%Currently requires that M is a KronProd consisting of 
%numeric or logical matrix class operators only.
%
%The constituent operands are combined into one large matrix using kron().
%Scalar operands are converted to a multiple of the identity matrix when
%doing so.




M=sparse(M);

M=full(M);
