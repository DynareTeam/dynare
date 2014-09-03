function M=ctranspose(M)
%Complex conjugate transpose of a KronProd class
%
% M' returns a KronProd object representing full(M)'

for ii=1:M.numops

 M.opset{ii}=M.opset{ii}' ;

end

tmp=M.rangesizes;
M.rangesizes=M.domainsizes;
M.domainsizes=tmp;
