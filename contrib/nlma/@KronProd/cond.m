function val=cond(M,varargin)
%COND method for KronProd. Syntax is the same as for  cond@double()

if anym(M.domainsizes~=M.rangesizes) & ~isempty(varargin)
  if varargin{1}~=2
   error('Domain and Range are of unequal dimensions. Use the 2 norm.')
  end
end

if M.scalarcoeff==0,
  val=inf;
  return;
end


Conds=zeros(1,M.numops);
for ii=1:M.numops

 Conds(ii)=cond(full(M.opset{ii}),varargin{:});
 
end


Conds=Conds(M.opinds);

if nargin>1 && isequal(varargin{1},'fro')
 Conds(M.scalarmask)=Conds(M.scalarmask).*(M.domainsizes(M.scalarmask));
end

val=prod(Conds);