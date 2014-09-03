function out=mldivide(L,R)
%MLDIVIDE method for KronProd.
%
%In the explanations below, it will be helpful to be familiar with the class
%method full(M), whcih  converts a KronProd object M to the full numeric matrix 
%that it represents. See "help KronProd/full" and "help KronProd/sparse"
%
%
%(1) K=c\M  where M is a KronProd object and c is a scalar will return a
%KronProd object K representing c\full(M).
%
%
%(2) M=M1\M2 where M1 and M2 are both KronProd objects whose operands M1.opset and
%M2.opset are compatibly sized will attempt to construct a KronProd M by
%applying the \ operator to corresponding operands of M1 and M2. E.g., if
%
%M1 represents (A kron B)
%M2 represents (C kron D)
%
%then M1\M2 tries to construct a KronProd representing  (A\C kron B\D). 
%If the operands are sized compatibly with this, full(M) is equivalent to
%full(M1)\full(M2).
%
%
%(3) X=M\Y where M is a KronProd object and Y is a logical or numeric matrix 
%is numerically equivalent to X=full(M)\Y when full(M) and Y are compatibly sized, i.e., 
%when size(Y,1)=size(full(M),1)=prod(M.rangesizes). Of course, the computation exploits the 
%Kronecker product structure and is typically much faster than X=full(M)\Y.
%
%Additionally, if the columns of Y are reshaped so that size(Y)=[M.rangesizes, N], 
%then X=M\Y produces the same result as before but with the columns of X
%accordingly reshaped so that size(X)=[M.domainsizes, N]. 
%In other words, M\Y attempts to invert the operation Y=M*X when M acts as 
%an operator mapping arrays of size M.domainsizes to arrays of size M.rangesizes.




  if isnum(L) %means K=c\M

     c=L; M=R;   
     
     out=M;
     out.scalarcoeff=full(double(out.scalarcoeff/c));
     
    return;
    
  elseif isa(L,'KronProd')  &  isa(R,'KronProd') %means M=M1\M2

       M1=L; M2=R;

       if ~isequal(M1.rangesizes,M2.rangesizes), error 'Incompatible sizes'; end
       
       M2_opset=M2.opset;

       for ii=find(M2.scalarmask)
           idx=M2.opinds(ii);
           M2_opset{idx}=M2_opset{idx}*eye(M2.domainsizes(ii));
       end
       
       if isequal(M1.opinds, M2.opinds)

        opinds=M2.opinds; zzz=unique(opinds);
        opset=cell(1, length(zzz));   
        for ii=zzz,  opset{ii}=M1.opset{ii}\M2_opset{ii};  end

       else

        opinds=1:length(M2.domainsizes);
        opset=cell(size(opinds));
        for ii=opinds,  
            opset{ii}=M1.opset{M1.opinds(ii)}\M2_opset{M2.opinds(ii)};  
        end

       end
       
        domainsizes=M2.domainsizes;
        scalarcoeff=M1.scalarcoeff\M2.scalarcoeff;
        
        out=KronProd(opset,opinds,domainsizes,scalarcoeff);
        
        return;	

  elseif ~(isa(L,'KronProd')  &  ( isnumeric(R)| islogical(R)) )
    error('Undefined mldivide operation involving KronProd');
  end

%%%%

M=L; X=R;

restoresingle=false;

if isa(X,'single')
 for ii=M.opinds

    if issparse(M.opset{ii})
     X=double(X);
     restoresingle=true;
     break;
    end

 end
end
%%%%



numel_domain=prod(M.domainsizes);
numel_range=prod(M.rangesizes);


[wascolumnized,layers]=procShape(M.rangesizes,size(X));



%%


nn=M.maxdim;

current_dims=M.rangesizes;
if isnum(current_dims), current_dims(2)=1; end


scalar=M.scalarcumprods(nn)*M.scalarcoeff;



if scalar==0%Computationally easy case
    
    X=X(:);
    X(numel_range)=0;
    X=X(1:numel_range);
    X(1:numel_range)=inf;
    
else

        if numel_domain<=numel_range && scalar~=1
 
           X=scalar\X;

        end


%%LOOP OVER DIMENSIONS


     for ii=1:nn


        X=reshape(X,current_dims(ii),[]);

          if ~M.scalarmask(ii)%current operator is non-trivial operator
    
            X=M.opset{M.opinds(ii)}\X;

          end



        X=X.';



     end


    %%
    if numel_domain>numel_range && scalar~=1
 
           X=scalar\X;

    end
        
end%scalar=0

if layers>1,
     X=reshape(X,layers,[]).';
end


if wascolumnized
  X=reshape(X,numel_domain,[]);   
else
  X=reshape(X,[M.domainsizes,layers]);
end

if restoresingle
  X=single(X);
end

out=X;

