function varargout=subsref(M,S)
%SUBSREF mthod for KronProd class
%
%M.fieldname   - accesses fieldname
%
%M{DIMS}       - Returns the list of operators applied to the directions specified in DIMS
%
%M(DIMS)       - constructs a new KronProd which is the restriction of
%                M to dimensions listed in DIMS


opsetflag=false; 

switch S(1).type

case '.'

  str=indexstr(S);

  varargout=eval(['{M' str  '};']);  

  return

case '{}'

 inds=[S(1).subs{1}];
 inds = M.opinds(inds);
 B=M.opset(inds);
 if isnum(inds), B=B{:}; end

case '()'

 new_domainsizes=enrow(M.domainsizes( S(1).subs{1}  )); 
 new_opinds=enrow(M.opinds( S(1).subs{1}  )); 
 if islogical(new_opinds), new_opinds=find(newopinds); end
 if ~isvec(new_opinds), error('M(...) can only take a vector'); end

 ww=unique(new_opinds);
 new_numops=length(ww);

 tt=cell(1,max(ww));
 yy=encell(1:new_numops);
 [tt{ww}]=deal(yy{:});

 new_opinds= [tt{new_opinds}];


 B=KronProd(M.opset(ww), new_opinds,new_domainsizes,M.scalarcoeff);

end


if length(S)>1,

  B=subsref(B,S(2:end));

end


if iscell(B),

 if length(S)==1 & opsetflag
   varargout{1}=B; 
 else
  varargout=B;
 end


else
 varargout{1}=B;
end



function str=indexstr(S)
%Given a subscripting structure array S, such as output by substruct() and
%used in indexing methods SUBSREF and SUBSASGN, a single string STR is 
%returned bearing a subscript  expression equivalent to that expressed 
%by S. This is useful when combined  with EVAL() to implement overloaded 
%SUBSASGN and SUBSREF methods based on 
%already-defined subscript operations.
%
%Usage:
%
% str=indexstr(S)
%
%EXAMPLE:
%
%         >>a.g={1 2 3};
%         >>M=substruct('.', 'g', '{}', {':'}); 
%         >>str=indexstr(M)
% 
%           str =
% 
%            .g{M(2).subs{:}}
%           
%          >>eval(['a' str])   %equivalent to a.g{:}
% 
%             ans =
% 
%                  1
% 
% 
%             ans =
% 
%                  2
% 
% 
%             ans =
% 
%                  3
%
%NOTE: When it exists, the INPUTNAME of S is used in generating STR. Otherwise,
%      this name defaults to 'S'
%
%
%See also SUBSTRUCT, INPUTNAME


name=inputname(1);

%DEFAULT
if isempty(name), name='S'; end


%%
str='';

for ii=1:length(S)

 switch S(ii).type

 case '.'



   str=[str '.' S(ii).subs];

 case '()'

   insert=[name '(' num2str(ii) ')'];

   str=[str '(' insert '.subs{:})'];

 case '{}'

   insert=[name '(' num2str(ii) ')'];

   str=[str '{' insert '.subs{:}}'];

 end

end

