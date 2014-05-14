function [varargout]=eig(M1,varargin);
%EIG method for KronProd class
%
%The eigendecomposition of a Kronecker product is distributive across the
%operands of the product, if the operands are square.
%
%E = EIG(M,options{:}) tries to exploit this by applying Ei=eig@double() to each
%operand of the KronProd object, M. Then E is a KronProd object representing the
%Kronecker product over the Ei.
%
%Similarly, [V,D] = EIG(M,options{:}) produces a KronProds V of eigenvector matrices
%and D of diagonal eigenvalue matrices.
%
%WARNING: Any of the above outputs can be converted from KronProd objects
%to numeric form using full() or sparse(), but the result will be the same as 
%output of eig(full(M)) only up to a re-ordering of the eigenvalues/vectors.


    
GenProb=false;
if length(varargin)>0 && ~ischar(varargin{1})
    
    GenProb=true;
    M2=varargin{1};
    
    if ~isa(varargin{1},'KronProd') || ~isequal(M1.domainsizes,M2.domainsizes)
     error 'Generalized eigenvalue problem only operable between two KronProds of equal-sized operands.'
    end
    
end





%%
nn=nargout;

if ~GenProb
 
        Dcoeff=M1.scalarcoeff;
        opinds=M1.opinds;
        scalarmask=M1.scalarmask;
        eyemask=M1.eyemask;
        domainset=M1.domainset;
        maxdim=M1.maxdim;
        domainsizes=M1.domainsizes;   
        
        N=M1.numops;
        
       ArgsOut=cell(nn,N);
    
        for ii=1:N %loop over opset

             [ArgsOut{:,ii}]=eig(M1.opset{ii},varargin{:});

        end


else%generalized eigenvalue problem
   
    Dcoeff=M1.scalarcoeff/M2.scalarcoeff;
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

        [ArgsOut{:,ii}]=eig(M1operand,M2operand,varargin{2:end});


      end

     
    
end
        
        
        
        
if nn>1
   
    %Domains=cellfun('size',ArgsOut,2);
    Domains=cellfun(@(x)size(x,2),ArgsOut);
    Domains=Domains(:,opinds);
    Domains(:,scalarmask)=repmat(domainsizes(scalarmask),nn,1);

    for jj=1:nn

        s=1; if ii==2, s=Dcoeff; end
        varargout{jj}=KronProd(ArgsOut(jj,:), opinds,Domains(jj,:),s);

    end
    
else%n==1
   
    for ii=1:N,
        
        if eyemask(ii)
           ArgsOut{ii}=repmat(ArgsOut{ii},domainset(ii),1); 
        end
    
    end 
    
    varargout{1}=KronProd(ArgsOut,opinds,ones(1,maxdim),Dcoeff);
    
end
