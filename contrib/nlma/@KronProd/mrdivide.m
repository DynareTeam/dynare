function out=mrdivide(L,R)
%MRDIVIDE for KronProd class
%
%
%In the explanations below, it will be helpful to be familiar with the class
%method full(M), whcih  converts a KronProd object M to the full numeric matrix 
%that it represents. See "help KronProd/full" and "help KronProd/sparse"
%
%
%(1) K=M/c  where M is a KronProd object and c is a scalar will return a
%KronProd object K representing full(M)/c.
%
%
%(2) M=M1/M2 where M1 and M2 are both KronProd objects whose operands M1.opset and
%M2.opset are compatibly sized will attempt to construct a KronProd M by
%applying the / operator to corresponding operands of M1 and M2. E.g., if
%
%M1 represents (A kron B)
%M2 represents (C kron D)
%
%then M1/M2 tries to construct  a KronProd representing  (A/C kron B/D). 
%If the operands are sized compatibly with this, full(M) is equivalent to
%full(M1)/full(M2).
%
%
%(3) Yt=Xt/M where M is a KronProd object and Yt is a logical or numeric matrix 
%is numerically equivalent to Yt=Xt/full(M) when full(M) and Xt are compatibly sized, i.e., 
%when size(Xt,2)=size(full(M),2)=prod(M.domainsizes). Of course, the computation exploits the 
%Kronecker product structure and is typically much faster than Yt=Xt/full(M).
%
%Additionally, if the rows of Xt are reshaped so that size(Xt)=[N,M.domainsizes], 
%then Yt=Xt/M produces the same result as before but with the rows of Yt
%accordingly reshaped so that size(Yt)=[N,M.rangesizes]. 
%In other words, Xt/M attempts to invert the operation Yt*M=Xt in a way
%that is compatible with the size restrictions described in "help KronProd/mtimes".


if isnum(R) %Means  K=M/c
    
    
    M=L; c=R;  %Solving 
    
    
    out=M*(1/c); return
    
elseif isa(L,'KronProd')  &  isa(R,'KronProd')  %Means M=M1/M2
    
    M1=L; M2=R; 
    
    out=((M2.')\(M1.')).';
    
    return;
    
    
elseif isnumeric(L)|islogical(L), %Means Y=X/M, i.e. solve Y*M=X
    
       X=L; M=R; 
       
    
      [wasrowshaped,layers]=procShape(M.domainsizes,size(X),1);
    
      Xt=reshape(X,layers,[]).';
      
      out=(M.'\Xt).';
      
      if ~wasrowshaped,
       sz=M.rangesizes;
       if layers>1, sz=[layers,sz]; end
       out=reshape(out,sz);
      end
      
      return
    
    
else
     error('Undefined usage of KronProd/mrdivide');
end


