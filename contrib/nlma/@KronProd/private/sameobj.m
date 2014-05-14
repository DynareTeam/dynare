function bool=sameobj(L,R)
%Checks to see if L and R are the same numerical object.
%
%To return true, L and R must be equal in size, class, and all their entries

bool=strcmp(class(L),class(R));
if ~bool, return; end

bool=allm(size(L)==size(R));
if ~bool, return; end

bool=allm(L==R);
