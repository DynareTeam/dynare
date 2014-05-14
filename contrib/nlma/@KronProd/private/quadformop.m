function M=quadformop(L,R)
%M is the matrix such that [M*v](:) = [L*diag(v)*R](:) 
%
%Currently, M is always expressed as sparse.
%
%If R is ommitted, it is assumed that R=L'

if nargin<2
 R=L';
end

nn=size(L,2);

%M=zeros( size(L,1)*size(R,2) , nn);
M=sparse( size(L,1)*size(R,2) , nn);

for ii=1:nn %THIS WOULD BE BETER AS MEX

 cc=L(:,ii) * R(ii,:); 
 M(:,ii)=cc(:);

end
