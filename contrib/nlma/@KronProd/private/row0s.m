function out=row0s(n,outtype)
% Outputs a row vector of N zeros


if nargin<2

 out=zeros(1,n);

else

 if strcmp(outtype,'logical')

  out=false(1,n);

 elseif vernum>=7

  out=zeros(1,n,outtype);

 end



end
