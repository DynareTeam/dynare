function M=inv(M)
%INV method for KronProd class
%
%  Y=inv(M)
%
%Each operand of M must be square, non-singular.

invopset=cell(1,M.numops);

for ii=1:M.numops

 invopset{ii}= inv(M.opset{ii});

end


opinds=M.opinds;
scalarcoeff=1/M.scalarcoeff;
domainsizes=M.rangesizes;

M=KronProd(invopset,opinds,domainsizes,scalarcoeff);