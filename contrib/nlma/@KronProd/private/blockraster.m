function Y=blockraster(X,rowdims,coldims)

if length(rowdims)~=length(coldims)
 error('ROWDIMS vector and COLDIMS must have equal length');
end

nn=dimprod(X);
Y=reshape(X,rowdims(1),[]); clear X;

if length(rowdims)<2
 if coldims~=size(Y,2); error 'Incosistent COLDIM'; end
 return;
end

rowdims=enrow(rowdims);  
coldims=enrow(coldims);

rowprods=cumprod(rowdims);
colprods=cumprod(coldims);

stripwidths=nn./rowprods;
numhorztiles=stripwidths./colprods;


for ii=2:length(rowdims)


  Y=mat2cell(Y,rowprods(ii-1),horzcp(colprods(ii-1),numhorztiles(ii-1)) );
  Y=reshape(Y,rowdims(ii),[]);
  Y=cell2mat(Y);


end
