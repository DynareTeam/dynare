function val=norm(M,varargin)
%NORM method for KronProd. Syntax is the same as for  cond@norm()


val=abs(M.scalarcoeff);

if val==0, return; end

Norms=zeros(1,M.numops);

for ii=1:M.numops

 Norms(ii)=norm(full(M.opset{ii}),varargin{:});

end

Norms=Norms(M.opinds);

if nargin>1 && isequal(varargin{1},'fro')
 Norms(M.scalarmask)=Norms(M.scalarmask).*sqrt(M.domainsizes(M.scalarmask));
end
 
 
val=val*prod(Norms);
