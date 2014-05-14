function out=mtimes(L,R)
%MTIMES method for KronProd. (Revised Dec. 4, 2009)
%
%In the explanations below, it will be helpful to be familiar with the class
%method full(M), which  converts a KronProd object M to the full numeric matrix 
%that it represents. See "help KronProd/full" and "help KronProd/sparse"
%
%
%(1) K=c*M or K=M*c where M is a KronProd object and c is a scalar will return a
%KronProd object K representing c*full(M).
%
%
%(2) M=M1*M2 where M1 and M2 are both KronProd objects whose operands M1.opset and
%M2.opset are compatibly sized will attempt to construct a KronProd M by
%multiplying corresponding operands of M1 and M2. E.g., if
%
%M1 represents (A kron B)
%M2 represents (C kron D)
%
%then M1*M2 tries to construct a KronProd representing (A*C kron B*D). 
%If the operands are sized compatibly with this, full(M) is equivalent to 
%full(M1)*full(M2).
%
%
%(3) Y=M*X where M is a KronProd object and X is a logical or numeric matrix 
%is numerically equivalent to Y=full(M)*X when full(M) and X are compatibly sized, i.e., 
%when size(X,1)=size(full(M),2)=prod(M.domainsizes). Of course, the computation exploits the 
%Kronecker product structure and is typically much faster than Y=full(M)*X.
%
%Additionally, if the columns of X are reshaped so that size(X)=[M.domainsizes, N], 
%then Y=M*X produces the same result as before but with the columns of Y
%accordingly reshaped so that size(Y)=[M.rangesizes, N]. In other words,
%M acts as an operator mapping arrays of size M.domainsizes to 
%arrays of size M.rangesizes.
%
%
%(4) X=Y*M where M is a KronProd object and Y is a logical or numeric matrix 
%is numerically (but again not algorithmically) equivalent to X=Y*full(M) when 
%full(M) and Y are compatibly sized, i.e., when size(Y,2)=size(full(M),1)=prod(M.rangesizes). 
%
%Additionally, if the rows of Y are reshaped so that size(Y)=[N,M.rangesizes], 
%then X=Y*M produces the same result as before but with the rows of X
%accordingly reshaped so that size(X)=[N,M.domainsizes].



  if isnum(L) %means K=c*M

     c=L; M=R;   
     
     out=M;
     out.scalarcoeff=full(double(out.scalarcoeff*c));
     
    return;

  elseif isnum(R) %means K=M*c

     M=L; c=R;
      
     out=M;
     out.scalarcoeff=full(double(out.scalarcoeff*c));
     return;

  elseif isa(L,'KronProd')  &  isa(R,'KronProd') %Means M=M1*M2;

   M1=L; M2=R;
      
   if ~isequal(M1.domainsizes,M2.rangesizes), error 'Incompatible sizes'; end
   
   
   domainsizes=M2.domainsizes;
   scalarcoeff=M1.scalarcoeff*M2.scalarcoeff;

   if isequal(M1.opinds, M2.opinds)

    opinds=M2.opinds; zzz=unique(opinds);
    opset=cell(1, length(zzz));   
    for ii=zzz,  opset{ii}=M1.opset{ii}*M2.opset{ii};  end
   
   else

    opinds=1:length(M2.domainsizes);
    opset=cell(size(opinds));
    for ii=opinds,  opset{ii}=M1.opset{M1.opinds(ii)}*M2.opset{M2.opinds(ii)};  end
     
   end

    out=KronProd(opset,opinds,domainsizes,scalarcoeff);
    
   return;	

  elseif isnumeric(L)|islogical(L)% means Y=X*M
      
      X=L; M=R;
      
      [wasrowshaped,layers]=procShape(M.rangesizes,size(X),1);
    
      Xt=reshape(X,layers,[]).';
      
      out=(M.'*Xt).';
      
      
      if ~wasrowshaped,
       sz=M.domainsizes;
       if layers>1, sz=[layers,sz]; end
       out=reshape(out,sz);
      end
      

      return
    
  elseif ~(isa(L,'KronProd')  &  (isnumeric(R)| islogical(R)) )
    error('Undefined mult. operation involving KronProd');
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

%wasvector=iscol(X);

numel_domain=prod(M.domainsizes);
numel_range=prod(M.rangesizes);



%%%%

[wascolumnized,layers]=procShape(M.domainsizes,size(X));




%%
nn=M.maxdim;

current_dims=M.domainsizes;
if isnum(current_dims), current_dims(2)=1; end


scalar=M.scalarcumprods(nn)*M.scalarcoeff;



if scalar==0%Computationally easy case
    
    X=X(:);
    X(numel_range)=0;
    X=X(1:numel_range);
    X(1:numel_range)=0;
    
else

        if numel_domain<=numel_range && scalar~=1
 
           X=scalar*X;

        end


%%LOOP OVER DIMENSIONS


     for ii=1:nn


        X=reshape(X,current_dims(ii),[]);

          if M.scalarmask(ii)%current operator represented by scalar

              X=X.'; %nothing else necessary
              
          elseif issparse(M.opset{M.opinds(ii)}) && ~issparse(X)

           %FULL*SPARSE MUCH FASTER THAN SPARSE*FULL

           %X=(X.' * M.opset{M.opinds(ii)}.').';
           
           X=(X.' * M.opset{M.opinds(ii)}.');  %Note: no tranpose necessary

          else

           X=M.opset{M.opinds(ii)} * X;
           X=X.';
           
          end

           %X=X.';

       



     end


    %%
    if numel_domain>numel_range && scalar~=1
 
           X=scalar*X;

    end
        
end%scalar=0


if layers>1,
     X=reshape(X,layers,[]).';
end


if wascolumnized
  X=reshape(X,numel_range,[]);   
else
  X=reshape(X,[M.rangesizes,layers]);
end


if restoresingle
  X=single(X);
end

out=X;

