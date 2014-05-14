function N=numel(M,varargin)
%NUMEL for a KronProd


if ~isempty(varargin)

 N=numel(M.opinds,varargin{:});

else

 N=1;

end



