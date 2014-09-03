function n=vernum
%Gets the version number up to one decimal place.

sss=version;


dec=strfind(sss,'.');
sss=sss(1:dec(2)-1);
n=str2num(sss);
