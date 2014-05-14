function M=transpose(M)
%Transposes a KronProd object
%
% M.' returns a KronProd object representing full(M).'

for ii=1:M.numops

 M.opset{ii}=M.opset{ii}.' ;

end


tmp=M.rangesizes;
M.rangesizes=M.domainsizes;
M.domainsizes=tmp;
