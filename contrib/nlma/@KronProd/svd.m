function [varargout]=svd(M,varargin);
%SVD method for KronProd class
%
%The singular value decomposition of a Kronecker product is distributive across the
%operands of the product.
%
%[U,S,V] =SVD(M,options{:}) effectively applies [Ui,Si,Vi]=svd@double() to
%all operands of the KronProd, M. U is a KronProd object representing the Kronecker
%product over Ui and analogously for S and V.
%
%S=SVD(M,options{:}) produces a KronProd of the singular value vectors.
%
%WARNING: Any of the above outputs can be converted from KronProd objects
%to numeric form using full() or sparse(), but the result will be the same as 
%output of svd(full(M)) only up to a re-ordering of the singular values/vectors.





nn=nargout;
N=M.numops;

%%
for ii=1:N

     [ArgsOut{ii,1:nn}]=svd(M.opset{ii} ,varargin{:});
 
end

if nn>1

    %Domains=cellfun('size',ArgsOut,2);
    Domains=cellfun(@(x)size(x,2),ArgsOut);
    Domains=Domains(M.opinds,:).';
    Domains(:,M.scalarmask)=repmat( M.domainsizes(M.scalarmask) ,nn,1);

    for ii=1:nn

        s=1; if ii==2, s=M.scalarcoeff; end
        varargout{ii}=KronProd(ArgsOut(:,ii),M.opinds,Domains(ii,:),s);

    end
    
else
   
    for ii=1:M.maxdim,
        
        if M.scalarmask(ii)
           idx=M.opinds(ii);
           ArgsOut{idx}=repmat(ArgsOut{idx},M.domainsizes(ii),1); 
        end
    
    end 
    
    varargout{1}=KronProd(ArgsOut,M.opinds,ones(1,M.maxdim),M.scalarcoeff);
    
end


