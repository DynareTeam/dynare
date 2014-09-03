function varargout=size(M,varargin);
%SIZE method for KronProd class.
%
%     varargout=size(M,varargin)
%
%is equivalent to 
%
%     varargout=size(full(M),varargin)
%
%where full(M) is the representation of M as a full numeric matrix 
%(see "help KronProd/full"). 
%
%In other, words [P,Q]=size(M) returns P=prod(M.rangesizes) and
%Q=prod(M.domainsizes).




sz=[prod(M.rangesizes), prod(M.domainsizes)];

varargout=parse_sizevec(sz,nargout,varargin{:});
