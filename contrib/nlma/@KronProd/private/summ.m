function out=summ(M)
%SUMM(M) is the sum over the entries of the array M
%If M is  a non-numeric array, but supports {} subscripting
%then sum(M(:)){1} is returned.

out=(sum(M(:)));

if isnumeric(out);

 out=full(out);

else

 out=full(out{1});

end
