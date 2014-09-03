function M=sum(M,varargin);
%SUM method for KronProd class

if nargin<2
    dim=find(size(M)~=1,1);
else
    dim=varargin{1};
    varargin(1)=[];
end


domainsizes=M.domainsizes;
if dim==2, domainsizes(:)=1; end 


for ii=1:M.numops

  M.opset{ii}=  sum(M.opset{ii}, dim, varargin{:});

  if M.eyemask(ii)
  
      switch dim
   
            case 1

              M.opset{ii}=M.opset{ii}*ones(1,M.domainset(ii));  


            case 2
                
              M.opset{ii}=M.opset{ii}*ones(M.domainset(ii),1);  
                
      end
   end
end


M=KronProd(M.opset,M.opinds,domainsizes,M.scalarcoeff);
