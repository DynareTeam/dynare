function M=abs(M)
%ABS() method for @KronProd
%

absopset=cell(1,M.numops);

for ii=1:M.numops

 absopset{ii}= abs(M.opset{ii});

end


opinds=M.opinds;
scalarcoeff=abs( M.scalarcoeff  );
domainsizes=M.domainsizes;

M=KronProd(absopset,opinds,domainsizes,scalarcoeff);