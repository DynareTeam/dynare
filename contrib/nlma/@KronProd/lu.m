function [varargout]=lu(M,varargin);
%LU method for KronProd class
%
%When the operands of a Kronecker product are square matrices, the LU decomposition
%is distributive across the operands of the product, i.e., LU factors 
%of the Kronecker product can be obtained as the Kronecker product of the
%LU factors.
%
%For a KronProd object, M, the syntax
%
%  [varargout{1:N)] = lu(M,varargin{:}) 
%
%tries to exploit this by applying
%
%    [ArgsOut{i,1:N)] = LU@double(M.opset{i},varargin{:}) 
%
%to each individual operand M.opset{i} of the KronProd object, M. 
%
%For every j=1...N, each set varargout_opset{:,j} then become Kronecker product operands
%to the j-th output argument of the method, which is likewise a KronProd object.
%
%NOTE 1: If the operands are not square, the code will still carry out the
%processing described above, but that resulting factorization may not be 
%an LU factorization in the strict sense.
%In particular, full(L) and/or full(U) may only be block triangular.
%However, the result may still be useful for manipulating LU
%factors tensorially.
%
%NOTE 2: If even one operand of the Kronecker product represented by M is
%sparse, then all  operands will be processed as sparse and the
%input/output behavior of the method will mimic that of a sparse matrix.
%
%NOTE 3: The 'vector' input arguments from lu@double(...) are not currently supported.
%        
%NOTE 4: The single argument syntax Y=lu@double(X) is not supported here, currently.


mm=nargin;
nn=nargout; 
 if nn<=1, 
     error 'Single output argument case not supported'
 end
N=M.numops;


    if any( cellfun(@(x) isequal(x,'vector'), varargin) )
        error 'The ''vector'' input argument is not supported in this class method'
    end

if ~isequal(M.domainsizes,M.rangesizes)
warning('KronProd:LU:OperandSizes',...
       'For nonsquare Kronecker product operands, the resulting LU factors may only be block triangular') 
end

SparseFlag=any( cellfun(@(x) issparse(x), M.opset(M.opinds)) );

if SparseFlag, 
    M.opset=cellfun(@(x) sparse(x), M.opset,'uniformoutput',false); 
end

%%
for ii=1:N

     [ArgsOut{ii,1:nn}]=lu(M.opset{ii},varargin{:});

end



    %Domains=cellfun('size',ArgsOut,2);
    Domains=cellfun(@(x)size(x,2),ArgsOut);
    
    Domains=Domains(M.opinds,:).';
    Domains(:,M.scalarmask)=repmat( M.domainsizes(M.scalarmask) ,nn,1);
    
   
    for ii=1:nn
        s=1;
    
        if ii==2, s=M.scalarcoeff; end
        varargout{ii}=KronProd(ArgsOut(:,ii),M.opinds,Domains(ii,:),s);
    end   
        
  