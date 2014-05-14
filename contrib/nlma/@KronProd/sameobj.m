function bool=sameobj(L,R)
%Determines equality of two KronProds

bool= L.scalarcoeff==R.scalarcoeff;

bool= bool & L.maxdim==R.maxdim;


if ~bool, return; end


for ii=1:L.maxdim

 bool= bool & L.opset{L.opinds(ii)}==R.opset{R.opinds(ii)};

 if ~bool, return; end

end
