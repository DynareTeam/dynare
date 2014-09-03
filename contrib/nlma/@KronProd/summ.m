function val=sum(M);
%SUMM method for KronProd class

M=sum(M,1);
M=sum(M,2);

val=full(M);
