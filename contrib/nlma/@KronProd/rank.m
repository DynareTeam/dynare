function r=rank(M,varargin);
%  RANK  method for KronProd class.
%
%     RANK(M) provides an estimate of the number of linearly
%     independent rows or columns of full(M) where M is a KronProd object.
%     RANK(M,tol) is the number of singular values of M
%     that are larger than tol.
%     RANK(M) uses the default tol = max(size(M)) * eps(norm(M)).


s = full(svd(M));
if nargin==1
   tol = max(size(M)) * eps(max(s));
end
r = sum(s > tol);
