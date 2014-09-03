function bool=equaldims(size1,size2)
%Detect equal-sized arrays
%
%  bool=equaldims(size1,size2)
%
%Vectors size1 and size2 list the dimensions of arrays. They are 
%considered equal (bool=true) if they are identical after padding them
%as necessary with ones to make them equal length.




        size1=size1(:); size2=size2(:); 
        
        ii=find(size1~=1,1,'last'); jj=find(size2~=1,1,'last');
    
    

        bool=~isempty(ii==jj) && (ii==jj)&& all(  size1(1:ii)==size2(1:jj)  );
  
