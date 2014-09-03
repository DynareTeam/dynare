function M=pinv(M,varargin)
%PINV method for KronProd class
%
%  Y=pinv(M,varargin)
%
%Works by applying pinv@double(Mi,varargin) to each M.opset{i}.
%
%Sparse operands will be pre-converted to full (since pinv is only
%applicable to full matrices).
%
%Code will fail if pinv() fails on any operand.

invopset=cell(1,M.numops);

for ii=1:M.numops

 A=M.opset{ii};
 
 if issparse(A), A=full(A); end
    
 invopset{ii}= pinv(A,varargin{:});

end


opinds=M.opinds;
scalarcoeff=1/M.scalarcoeff;
domainsizes=M.rangesizes;

M=KronProd(invopset,opinds,domainsizes,scalarcoeff);