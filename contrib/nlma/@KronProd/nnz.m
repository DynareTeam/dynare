function val=nnz(M);
%NNZ method for KronProd class


if ~M.scalarcoeff
 val=0; return;
end

%%%

aa=row1s(M.numops);

for ii=1:M.numops

 aa(ii)=  nnz(M.opset{ii});

 if ~aa(ii), val=0; return; end

end


aa=aa(M.opinds);

aa(M.scalarmask)=M.domainsizes(M.scalarmask);

val=prod(aa);
