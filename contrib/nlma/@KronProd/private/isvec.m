function bool=isvec(a)
%tests if argument is a  vector. 
%An empty matrix  is considered a vector.

sz=size(a);
bool = ( nnz( sz>1  )  <= 1 ) & allm(sz>0) ;
bool=bool & allm(sz>0);
bool= bool | allm(sz==0);