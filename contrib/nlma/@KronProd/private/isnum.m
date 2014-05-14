function bool=isnum(a)
%Tests if argument is a numeric 1x1 object
%
%bool=isnum(a)
%
%Note: isnum(true)=true, isnum(NaN)=true

bool=is1x1(a) & (isnumeric(a) | islogical(a));

