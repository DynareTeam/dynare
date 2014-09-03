function Y=horzcp(X,nn)
%
%Horizontal copy.
%The output Y is the vertical concatenation of NN copies of 
%the input matrix X.


if nn<1; Y=[];return;end

if iscol(X)

  Y=X*row1s(nn,class(X));

else %the following seems to be twice as fast as using kron()

 [t,s]=size(X);
 q=nn*s;

 Y=zeros(t,q);
 aa=reshape(1:q ,[s,nn]);

 for ii=1:nn

 Y(:,aa(:,ii))=X; 

 end

end
