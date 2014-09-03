function B=loadobj(A)
%LOADOBJ method for KronProd

B=KronProd(A.opset,A.opinds,A.domainsizes,A.scalarcoeff);

% if isstruct(A) | issparse(A.scalarcoeff); 
%    B=KronProd(A.opset,A.opinds,A.domainsizes,A.scalarcoeff)  
% else
%     B=A; %No conversion necessary
% end