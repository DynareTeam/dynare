function M=kron(varargin)
%KRON method for @KronProd class
%
%USAGE:
%
%   K=kron(A,B,C,D,...)
%
%where A,B,C,D,... is a combination of KronProd objects and other object
%types. The result K will be a KronProd object representing their
%Kronecker product.
% 
%EXAMPLE:
%
%     Q=fft(eye(4));       %fft matrix
%
%     Qo=KronProd({Q},1);  %wrap Q in a KronProd object
% 
%     E=kron(Qo,Q);      %equivalent to both kron(Q,Q) and kron(Qo,Qo)
%
%
%     >>whos E
%       Name       Size            Bytes  Class       Attributes
% 
%       E         16x16              796  KronProd        
% 
%
%
%     >>E*A-fft2(A)  %compare numeric fft to implementation using E
% 
%     ans =
% 
%       1.0e-015 *
% 
%             0                  0             0.4441                  0          
%             0                  0                  0                  0          
%             0            -0.3331 - 0.1110i   0.1110            -0.3331 + 0.1110i
%             0                  0                  0                  0          


varargin=fliplr(varargin);


idxkp=cellfun('isclass',varargin,'KronProd');

scalarcoeff=prod(  cellfun(@(ob)ob.scalarcoeff, varargin(idxkp)  )  );

domainsizes=varargin;
domainsizes(idxkp)=cellfun(@(ob)ob.domainsizes, varargin(idxkp),'uniformoutput',0);
domainsizes(~idxkp)=cellfun(@(A)size(A,2), varargin(~idxkp),'uniformoutput',0);


opset=varargin;
opset(idxkp)=cellfun(@(ob)ob.opset, varargin(idxkp),'uniformoutput',0);
opset(~idxkp)=cellfun(@(A){A}, varargin(~idxkp),'uniformoutput',0);

Lengths=cellfun('length',opset);
Lengths=num2cell( cumsum([0, Lengths(1:end-1)]));

opinds=varargin;
opinds(idxkp)=cellfun(@(ob)ob.opinds, varargin(idxkp),'uniformoutput',0);
opinds(~idxkp)=cellfun(@(A) 1, varargin(~idxkp),'uniformoutput',0);


opinds=cellfun(@(a,b) a+b, opinds, Lengths, 'uniformoutput',0);

opset=[opset{:}];
opinds=[opinds{:}];
domainsizes=[domainsizes{:}];

M=KronProd(opset,opinds,domainsizes,scalarcoeff);