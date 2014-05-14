function out=issquare(a)
%tests if argument is a square array (of whatever type).


if (size(a,1)==size(a,2)) & (length(size(a))<3)

 out=1;

else

 out=0;

end
