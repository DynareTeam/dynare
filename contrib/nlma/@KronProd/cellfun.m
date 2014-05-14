function out=cellfun(hfunc, varargin)
%CELLFUN - overloaded cellfun method for KronProd class
%
%
%SYNTAX 1: 
%
%   out=cellfun(fun,K,'param1', value1,...)
%
%where fun is a string and K is a KronProd object is equivalent to
%
%   out=cellfun(hfunc,K.opset(K.opinds),'param1', value1,...)
%
%
%
%SYNTAX 2:
%
%   M=cellfun(fun,K)
%
%where fun is a function handle and K is a KronProd object, returns another
%KronProd object M such that M.opset{i} and M.scalarcoeff are obtained by 
%applying fun() to  K.opset{i) and K.scalarcoeff respectively.
%
%If any K.opset{i} is a scalar c representing c*speye(N), then K.opset{i}
%will be replaced by default with c*speye(N) before fun() is applied.
%
%
%
%SYNTAX 3:
%
%   M=cellfun(fun,K,'ExpandScalars',LogicalValue)
%
%is the same as syntax 2 when LogicalValue=1 (true). 
%
%When LogicalValue=0 (false), then if any K.opset{i} is a scalar c representing
%c*speye(N), the operand K.opset{i} will ***NOT*** be replaced with c*speye(N) when 
%fun() is applied. Instead, the code will assume that fun(c*speye(N)) =fun(c)*speye(N) 
%and the returned M will have M.opset{i}=fun(c). When this assumption is
%applicable, it obviously saves computation and memory to use this option.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Many operations on a Kronecker product are distributable across its operands, 
%a fact which many KronProd class methods exploit.
%
%For operations not already implemented as class methods, however,
%the KronProd cellfun() method gives a general way to apply operations
%distributively across the operands via syntaxes 2 and 3. Of course, it is still 
%permissible, and sometimes useful, to apply the cellfun() method to operations
%that are not truly distributable with respect to Kronecker multiplication.
%
%EXAMPLES:
%
%Consider the data
%
%  K=KronProd({[1 2; 3 4],pi},[1 2],[nan,2],4); 
%
%  Kfull=full(K),
% 
%     Kfull =
% 
%        12.5664   25.1327         0         0
%        37.6991   50.2655         0         0
%              0         0   12.5664   25.1327
%              0         0   37.6991   50.2655
% 
%
%
%Example 1: The function
%
%   fun= @(A) diag(diag(A));
%            
%is distributable with respect to Kronecker multiplication when the operands
%of the Kronecker product are all square matrices. The function clearly also 
%satisfies the assumption fun(c*speye(N))=fun(c)*speye(N). 
%
%It is therefore appropriate to apply cellfun() with ExpandScalars=false,
%although the same output below would be obtained either way,
%
%         >> full( cellfun(fun,K,'ExpandScalars',0) )
% 
%         ans =
% 
%            12.5664         0         0         0
%                  0   50.2655         0         0
%                  0         0   12.5664         0
%                  0         0         0   50.2655
%
%
%This clearly agrees with fun(Kfull), as one can see from above.
%
%
%Example 2: Conversely, the function
%
%    fun= @(A) sum(A(:));
%            
%does not satisfy the assumption fun(c*speye(N))=fun(c)*speye(N), although
%it is distributable with respect to Kronecker multiplication. Thus, the
%correct result is only obtained when ExpandScalars=true,
%
%         fun(Kfull),
% 
%         ans =
% 
%           251.3274
% 
%         full( cellfun(fun,K,'ExpandScalars',1) )  %correct result
% 
%         ans =
% 
%           251.3274
% 
%         full( cellfun(fun,K,'ExpandScalars',0) ) %incorrect result
% 
%         ans =
% 
%           125.6637         0
%                  0  125.6637



argMap=cellfun('isclass',varargin,'KronProd');

nn=find(argMap,1,'last');

dataMap=varargin(argMap);
varargin=varargin(~argMap);



if nn~=1 | ~all(argMap(1:nn))

   error 'Currently 1 and only 1 KronProd object must be an argument to KronProd.cellfun (and it must be the 2nd argument)' 
   %In future releases, this restriction might be relaxed.


else    

     K=dataMap{1};
    
     
end

if ischar(hfunc) %SYNTAX 1

     out=cellfun(hfunc,K.opset(K.opinds),varargin{:});
     
elseif isa(hfunc,'function_handle') %SYNTAX 2 and 3

       
    
        try
         varargin=reshape(varargin,[],2);
        catch
           error 'Options must occur in pairs' 
        end


        idx=strmatch('uni', lower( varargin(:,1) ));
        varargin(idx,:)=[];
        idx=strmatch('error', lower( varargin(:,1) ));
        varargin(idx,:)=[];


   
         idx=strmatch('exp', lower( varargin(:,1) ));
         ExpandScalars=true; %default
         if ~isempty(idx) 
           ExpandScalars=varargin{idx,2};
           varargin(idx,:)=[];
         end

         OpsetIn=K.opset;
         if ExpandScalars

             for ii=1:length(K.eyemask)

                 if ~K.eyemask(ii), continue; end
                 OpsetIn{ii}=OpsetIn{ii}*speye(K.domainset(ii));

             end


         end





         OpsetOut=cellfun(hfunc,OpsetIn,varargin{:},'UniformOutput',0);

         %domainsizes=cellfun('size', OpsetOut,2);
         domainsizes=cellfun(@(x)size(x,2), OpsetOut);
         domainsizes=domainsizes(K.opinds);
         
         if ~ExpandScalars
            domainsizes(K.scalarmask)=K.domainsizes(K.scalarmask); 
         end

         M=KronProd(OpsetOut,K.opinds,domainsizes,hfunc(K.scalarcoeff));
         out=M;



 else

     error 'Argument "fun" expected to be function handle or string.'

 end
