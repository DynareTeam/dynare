function [columnized,layers]=procShape(TrueSize,CatSize,transposed)
%PROCSHAPE - A helper function for the KronProd class used to determine
%whether a set of arrays (to be tensorially transformed) were columnized and
%concatenated into a 2D matrix, or else concatenated in their
%true n-dimensional shape.
%
%  columnized=procShape(TrueSize,CatSize)
%
%in:
%
% TrueSize: the dimensions of a single operand of the transform when 
%           represented n-dimensionally
%
% CatSize: the dimensions of the concatenated array.
%
% transposed: Boolean flag indicating that the transform targets have been
%             concatenated in a tranpose fashion. Default=false.
%
%out:
%
% columnized: A boolean output. TRUE if the operands were columnized.
%     layers: number of concatenated transform targets

if nargin<3, transposed=0; end


p=prod(TrueSize);
layers=prod(CatSize)/p;


if ~transposed
    
    z=strfind(CatSize,TrueSize);
    
    if ~isempty(z) && z(1)==1;

     
        columnized=false;

    elseif CatSize(1)==p;

        
        columnized=true;

    else

        error 'Cannot determine reshaping rule'

    end
    
else%Transposed reshaping
   
  
    
    if CatSize(1)==layers && CatSize(2)==p && length(size(CatSize))==2
   
        columnized=true;
        
    elseif equaldims(CatSize, TrueSize), %implies 1 layer
        
        columnized=false;
        
    elseif equaldims( CatSize(2:end), TrueSize)

        columnized=false;

    else

        error 'Cannot determine reshaping rule'

    end
    
end
