function X=rdivide(M,X)
%RDIVIDE method for KronProd. 
%
% Y=rdivide(M,X) or Y=M./X
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
%then M./X tries to construct a KronProd representing (A./C kron B./D). If
%the operands are sized compatibly with this, full(Y) is equivalent to
%full(M)./full(X).


   if isnum(M) %implies X is KronProd
  
     scalar=M;
    
     for ii=1:X.numops, X.opset{ii}=1./X.opset{ii}; end
     
     X.scalarcoeff=full(double(scalar/X.scalarcoeff));
     
     return;

  elseif  isnum(X) %implies M is KronProd

     tmp=X;
     X=M; 
     X.scalarcoeff=full(double(X.scalarcoeff./tmp));
     return;

  elseif isa(M,'KronProd')  &  isa(X,'KronProd')

 
   
   try, 
       compat=all(X.domainsizes==M.domainsizes & X.rangesizes==M.rangesizes);
   catch
       compat=false;
   end
   if ~compat
    error 'A.*B where A and B are both KronProds requires that A and B have operands of matching sizes'
   end
       
   domainsizes=X.domainsizes;
   scalarcoeff=M.scalarcoeff./X.scalarcoeff;

   
   
      Mopset=M.opset;
   Xopset=X.opset;
   
   for ii=find(X.scalarmask)
    idx=X.opinds(ii);   
    Xopset{ii}=Xopset{ii}*eye(X.domainsizes(ii));
   end
   
    for ii=find(M.scalarmask)
    idx=M.opinds(ii);   
    Mopset{ii}=Mopset{ii}*eye(M.domainsizes(ii));
   end
   
   if isequal(M.opinds, X.opinds)

    opinds=X.opinds; zzz=unique(opinds);
    opset=cell(1, length(zzz));   
    for ii=zzz,  opset{ii}=Mopset{ii}./Xopset{ii};  end
   
   elseif X.maxdim==M.maxdim

    opinds=1:length(X.maxdim);
    opset=cell(size(opinds));
    for ii=opinds,  opset{ii}=Mopset{M.opinds(ii)}./Xopset{X.opinds(ii)};  end
     
   end

    X=KronProd(opset,opinds,domainsizes,scalarcoeff);
    return;	

  else
    error('Undefined rdivide operation involving KronProd');
  end

