function H=quadform(L,v,R)
%Evaluate a weighted Gram matrix based on @KronProds
%
% H=quadform(L,v,R)
%
%Constructs H= L * diag(v)* R, and returns H as a sparse matrix.
%
%Assumes L and R are BOTH KronProds (and further R defaults to L'
%if omitted) with corresponding Kronecker constituents that are compatible 
%for matrix multiplication. Only in this case do we know an efficient 
%construction method for H.


if nargin<3,
 R=L';
end

nn=L.maxdim;

M=rowc(nn);


for ii=1:nn


 lblock=L.opset{L.opinds(ii)};
 rblock=R.opset{R.opinds(ii)};

 if isnum(lblock), lblock=speye(L.domainsizes(ii))*lblock; end
 if isnum(rblock), rblock=speye(R.domainsizes(ii))*rblock; end

 M{ii}=quadformop(lblock,rblock);

end


M=KronProd(M,1:nn,L.domainsizes); %THE 'KRONECKER QUADFORMOP'

%%

scalar=L.scalarcoeff*R.scalarcoeff;

v=v(:); 

if scalar~=1, v=scalar*v;  end %Jan. 09 -- Forgot about scalar coeffs!!

if ~issparse(v), v=sparse(double(v)); end

H=M*v;

H=blockraster(H, L.rangesizes,R.domainsizes);

if ~issparse(H), H= sparse(H); end
