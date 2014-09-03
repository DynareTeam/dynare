function out=row1s(n,outtype)
% Outputs a row vector of N ones


if nargin<2

 out=ones(1,n);

else

 if strcmp(outtype,'logical')

  out=true(1,n);

 elseif vernum>=7

  out=ones(1,n,outtype);

 end



end
