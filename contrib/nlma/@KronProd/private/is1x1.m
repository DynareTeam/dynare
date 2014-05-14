function out=is1x1(a)
%tests if argument is a 1 x 1 object


if (size(a,1)==1) & issquare(a)

 out=1;

else

 out=0;

end

out=logical(out);
