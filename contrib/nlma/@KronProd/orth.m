function Q=orth(M)
%ORTH method for the KronProd class
%
% Q = ORTH(M) produces a KronProd object satisfying Q'*Q = I and the columnspace of full(Q) 
% is the same as the column space of full(M). 
%
% The ORTH decomposition of a Kronecker product is distributive across the operands of the 
% product. This class method tries to exploit this by applying 
% orth@double() to each member of M.opset.
%
% Note: The number of columns of full(Q) may not agree with rank(Q). This is because
% the rank and column space of each M.opset{i} is individually analyzed here, 
% unlike rank(Q) which uses the svd of the entire object M
% to measure the rank. The latter is done to be consistent with
% rank@double().





for ii=1:M.numops
       
     opset{ii}=orth(M.opset{ii});
      
end

domainsizes=M.domainsizes;

for jj=1:M.maxdim
    if ~M.scalarmask(jj)
      domainsizes(jj)=size(opset{jj},2);
    end
end

Q=KronProd(opset,M.opinds,domainsizes,+1);