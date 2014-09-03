function [varargout]=qr(M1,varargin);
%QR method for KronProd class
%
%When the operand matrices of a Kronecker product have more columns than 
%rows, the QR decomposition is distributive across the operands of 
%the product, i.e., QR factors of the Kronecker product can be obtained as 
%the Kronecker product of the QR factors.
%
%For a KronProd object, M, the syntax
%
%  [varargout{1:N)] = qr(M,varargin{:}) 
%
%tries to exploit this by applying
%
%    [ArgsOut{i,1:N)] = qr@double(M.opset{i},varargin{:}) 
%
%to each individual operand M.opset{i} of the KronProd object, M. 
%
%For every j=1...N, each set varargout_opset{:,j} then become Kronecker product operands
%to the j-th output argument of the method, which is likewise a KronProd object.
%
%
%NOTE 1: If some operands have more columns than rows, the code will still carry out the
%processing described above, but that resulting factorization may not be 
%a QR factorization in the strict sense.
%In particular, full(R) may only be block upper triangular.
%However, Q will be unitary and the resulting factorization may still be useful if 
%one wishes to manipulate QR factors tensorially.
%
%NOTE 2: Unlike the method qr@double(...), the economy mode syntax [Arg1,Arg2,E]=QR(...,0) 
%always produces a permutation KronProd E, rather than a vector of permutation indices. 
%However, the operands of E will be sparse or, where possible, scalars standing in for
%identity matrices.
%        
%NOTE 3: If even one operand of the Kronecker product represented by M is
%sparse, then all  operands will be processed as sparse and the
%input/output behavior of the method will mimic that of a sparse matrix.
%
%NOTE 4: The single argument syntax Y=qr@double(X) is never distributable
%across Kronecker products when X is a full matrix. Hence, it is only
%supported currently when the KronProd object M is being processed as sparse
%(as outlined in NOTE 3).



mm=nargin;
nn=nargout; 


SparseFlag=any( cellfun(@(x) issparse(x), M1.opset(M1.opinds)) );

 if nn<=1 & ~SparseFlag, 
     error 'Single output argument case not supported for type full'
 end


if SparseFlag, 
    M1.opset=cellfun(@(x) sparse(x), M1.opset,'uniformoutput',false);  
end


if any(M1.domainsizes<M1.rangesizes)
   warning('KronProd:QR:OperandSizes',...
           'The QR factorization may have only block upper triangular R if some operands have more rows than columns.') 
end



EconomyFlag=any( cellfun(@(x) isequal(x,0), varargin) ) & (nn>=3);

QR_AB_FLAG=false;  M2isKronProd=false;
if length(varargin)>0 && ~isequal(varargin{1},0) && SparseFlag
    
    QR_AB_FLAG=true;
    
    if nn<2, error 'Must have at least 2 output args'; end
    
    M2=varargin{1};
    M2isKronProd=isa(varargin{1},'KronProd');
    
    if M2isKronProd && ~isequal(M1.rangesizes,M2.rangesizes)
     error 'QR(M1,M2,...) syntax between two KronProds M1 and M2 requires correponding operands to have equal number of rows.'
    end
    
    if ~M2isKronProd
       B=varargin{1}; 
       varargin(1)=[]; 
    end
    

end







%%

if ~(QR_AB_FLAG && M2isKronProd)
    

        opinds=M1.opinds;
        scalarmask=M1.scalarmask;
        eyemask=M1.eyemask;
        domainset=M1.domainset;
        maxdim=M1.maxdim;
        domainsizes=M1.domainsizes;   
        
        N=M1.numops;
    
        ArgsOut=cell(nn,N);
        
        
        for ii=1:N %loop over opset

           [ArgsOut{:,ii}]=qr(M1.opset{ii},varargin{:});

           if EconomyFlag & ~eyemask(ii)
              
               ee=ArgsOut{3,ii};
               E=speye(length(ee));
               E=E(:,ee);
                  if SparseFlag, E=E.'; end
               ArgsOut{3,ii}=E;
               
           end
             
             
        end

        
       
        
else%M2isKronProd
  

    scalarmask=M1.scalarmask & M2.scalarmask;
    maxdim=M1.maxdim;
    domainsizes=M1.domainsizes;   
    
    if isequal(M1.opinds, M2.opinds) %Can loop over minimized opset

           opinds=M1.opinds; 
           zzz=unique(opinds); 
 
             
           M1opset=M1.opset(zzz); 
           M1eyemask=M1.eyemask(zzz);
           M1domainset=M1.domainset(zzz);

           M2opset=M2.opset(zzz); 
           M2eyemask=M2.eyemask(zzz);
           M2domainset=M2.domainset(zzz);           
             
             
    else
        
          opinds=1:length(M2.domainsizes);
          M1opinds=M1.opinds;
          
          
           M1opset=M1.opset(M1.opinds); 
           M1eyemask=M1.eyemask(M1.opinds);
           M1domainset=M1.domainset(M1.opinds);

           M2opset=M2.opset(M2.opinds); 
           M2eyemask=M2.eyemask(M2.opinds);
           M2domainset=M2.domainset(M2.opinds);           
             
        
    end 
    

          
       eyemask=M1eyemask & M2eyemask;
       domainset=M1domainset;
       
       N=length(M1opset);
             
       ArgsOut=cell(nn,N);
      
      for ii=1:N,  

        M1operand=M1opset{ii};
          if ~eyemask(ii) & M1eyemask(ii)
             M1operand=M1operand*eye(M1domainset(ii));
          end

        M2operand=M2opset{ii};
          if ~eyemask(ii) & M2eyemask(ii)
             M2operand=M2operand*eye(M2domainset(ii));
          end

        zflag=false; 
        if isequal(M2operand,0) %could be confused with economy flag
         zflag=true;  
         M2operand=1;
        end
          
        [ArgsOut{:,ii}]=qr(M1operand,M2operand,varargin{2:end});

        if zflag, ArgsOut{1,ii}=0; end

           if EconomyFlag & ~eyemask(ii)
              
               ee=ArgsOut{3,ii};
               E=speye(length(ee));
               E=E(:,ee);
               if SparseFlag, E=E.'; end
               ArgsOut{3,ii}=E;
               
           end
 
      end
   
 
end

     


%%Post process operands

   
    %Domains=cellfun('size',ArgsOut,2);
    Domains=cellfun(@(x)size(x,2),ArgsOut);
    
    Domains=Domains(:,opinds);
    Domains(:,scalarmask)=repmat(domainsizes(scalarmask),nn,1);

    for jj=1:nn

        s=1; 
        if jj==2 |(jj==1 & nn==1 & SparseFlag), %the R scalarcoeff
            s=M1.scalarcoeff;
        elseif jj==1 && M2isKronProd  %the C scalar coeff
            s=M2.scalarcoeff;
        end,
        varargout{jj}=KronProd(ArgsOut(jj,:), opinds,Domains(jj,:),s);

    end
    
 
    
    if QR_AB_FLAG & ~M2isKronProd % and nn>1
        Q=varargout{1};
        C=Q'*B;
        varargout{1}=C;
    end
    
    

  