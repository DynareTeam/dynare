function [newobjset,newobjinds]=trimobjset(objset,objinds)
%A tool for compressing redundancies out of object sets.
%
%OBJSET is a list of objects, all of which are comparable via a SAMEOBJ()
%method, while OBJINDS is a vector of indices into this list.
%
%The output [NEWOBJSET,NEWOBJINDS] are equivalent, but no two objects in
%in NEWOBJSET are identical according to SAMEOBJ.

if ~isvar('objinds'), objinds=1:length(objset);end

if length(objset)<2, 
 newobjset=objset;
 newobjinds=objinds;
 return; 
end


objset=objset(objinds);
objinds=1:length(objset);

newobjset={objset{1}};
newobjinds=[1];

for ii=2:length(objinds)

  obj=objset{ii};

  for jj=1:length(newobjset)
  
   if sameobj(obj,newobjset{jj})
   
     newobjinds=[newobjinds,jj];
     break;
   end

   if jj==length(newobjset) %NON-REDUNDANT OBJECT
    newobjset=[newobjset,{obj}];
    newobjinds=[newobjinds,jj+1];
   end

  end

end
