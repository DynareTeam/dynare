function X=times(M,X)
%TIMES method for KronProd. 
%
% Y=times(M,X) or Y=M.*X
%
%At least one of M or X is a KronProd object. The other may be a scalar in which 
%case Y is an appropriately scaled KronProd object.
%     
%If M and X are both KronProd objects, then the operation is attempted on corresponding
%operands of M.opset and X.opset, and will only work if the sizes of the operands
%are compatible. E.g., if
%
%M1 represents (A kron B)
%M2 represents (C kron D)
%
%then M.*X tries to construct a KronProd representing  (A.*C kron B.*D). If
%the operands are sized compatibly with this, full(Y) is equivalent to full(M).*full(X).


  if isnum(M) %implies X is KronProd

    X.scalarcoeff=full(double(X.scalarcoeff.*M)); 
    return;

  elseif isnum(X) %implies M is KronProd

     tmp=X;
     X=M; 
     X.scalarcoeff=full(double(X.scalarcoeff.*tmp));
     return;

  elseif isa(M,'KronProd')  &  isa(X,'KronProd')

 
   
   try, 
       compat=isequal(X.domainsizes,M.domainsizes) & isequal(X.rangesizes, M.rangesizes);
   catch
       compat=false;
   end
   if ~compat
    error 'A.*B where A and B are both KronProds requires that A and B have operands of matching sizes'
   end
       
   domainsizes=X.domainsizes;
   scalarcoeff=M.scalarcoeff*X.scalarcoeff;

   if isequal(M.opinds,X.opinds)

    opinds=X.opinds; zzz=unique(opinds);
    opset=cell(1, length(zzz));   
    for ii=zzz,  opset{ii}=M.opset{ii}.*X.opset{ii};  end

   elseif M.maxdim==X.maxdim

    opinds=1:X.maxdim;
    opset=cell(size(opinds));
    for ii=opinds,  opset{ii}=M.opset{M.opinds(ii)}.*X.opset{X.opinds(ii)};  end
     
   end

   eyeproc=opinds(M.scalarmask|X.scalarmask);
   
   for ii=eyeproc
      opset{ii}=diag(diag(opset{ii})); 
   end
   
   
    X=KronProd(opset,opinds,domainsizes,scalarcoeff);
    return;	

  else
    error('Undefined times operation involving KronProd');
  end

  
