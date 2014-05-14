classdef KronProd
%KronProd class constructor.
%
%by Matt Jacobson, 
%Copyright, The University of Michigan, 2005
%
%
%This is a class for efficiently representing and manipulating
%N-fold Kronecker products of matrices (or of objects that behave like 
%matrices) in terms of their operands only. 
%
%Given matrices {A,B,C,D,...} and a scalar s, an object M of this class
%can be used to represent
%
%   Matrix = s * A kron B kron C kron D kron ...    (Eq. 1)
%
%where "A kron B" denotes kron(A,B), the Kronecker product of A and B. 
%Internally, however, M stores the data {s,A,B,C,D,...} separately,
%which  is typically far more byte-compact than numerically expanding out 
%the RHS of Eq. 1. 
%
%Furthermore, many mathematical manipulations of Kronecker products are more 
%efficient when done in terms of {s,A,B,C,D,...} separately than when done with 
%the explicit numerical form of M as a matrix. The class overloads a number 
%of methods and math operators in a way that exploits the Kronecker product
%structure accordingly.
%
%Among these methods/operators are: mtimes (*), times (.*) , tranpose (.') , 
%ctranpose (') ,  rdivide (./), ldivide (.\), mldivide (\), mrdivide (/), 
%inv, pinv, power, mpower, norm,  sum, cond,  eig, svd, abs, nnz, orth, chol, 
%lu, qr, cellfun, full, sparse,  ... 
%
%Some restrictions apply to these overloads. In particular, bi-operand math 
%operations  involving two KronProd objects, e.g. M1*M2, typically require 
%the operands of each KronProd to be of compatible sizes. However, I find these
%restrictions to be satisfied often in applications.
%Consult "help KronProd/methodname" for more info on each method.
%Optionally also, read krontest.m for demonstrations of their use.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXAMPLE #1:
%
%A primary application of this class is to efficiently perform separable
%tensorial operations, i.e., where a linear transform is applied to all
%columns of an array, then all rows, and so on.
%
%The following example is of a separable transformation of a 3D array X 
%that transforms all of its columns via multiplication with a non-square matrix A, then 
%transforms all rows by multiplication with B, then finally transforms all
%3rd-dimensional axes by multiplication with C.  Two approaches to this are compared.
%The first approach uses kron(). The second uses the KronProd class. Other operations
%are also shown for illustration purposes.
%
%Notice the orders of magnitude reduction both in CPU time and in memory
%consumption, using the KronProd object.
%
% 
%     %DATA
%     m=25; n=15; p=40;
%     mm=16; nn=n;  pp=10;
%     A=rand(mm,m); B=pi*eye(n); C=rand(pp,p); 
%     s=4; % a scalar
%     X=rand(m,n,p);
% 
% 
% 
% 
%     %METHOD I: based on kron()
%     tic;
% 
%      Matrix = s*kron(C,kron(B,A));
% 
%      y1 = Matrix*X(:);              %The tensorial transformation
%          y1=reshape(y1,[mm,nn,pp]); 
% 
%      z1 = Matrix.'*y1(:);
% 
%      w1 = Matrix.'\z1;
% 
% 
% 
%     toc;
%     %Elapsed time is 78.729007 seconds.
% 
% 
% 
% 
% 
%     %METHOD II: based on KronProd object
%     tic;
% 
%      Object = KronProd({A,pi,C},[1 2 3],[m,n,p],s); 
%         %
%         %equivalent to s*kron(C,kron(B,A)) - see USAGE section of help
% 
%        y2 =    Object*X;  
%                 % This operation could also have been implemented 
%                 % as y2=reshape( Object*X(:) , [mm,nn,pp]); 
%                      
% 
%                           
%        z2 = Object.'*y1;
% 
%        w2 = Object.'\z1;
% 
% 
%     toc
%     % Elapsed time is 0.003958 seconds.
% 
% 
% 
% 
%     %%%ERROR ANALYSIS
%     PercentError=@(x,y) norm(x(:)-y(:),2)/norm(x(:),'inf')*100;
% 
%         PercentError(y1,y2), % = 3.0393e-012
% 
%         PercentError(size(y1),size(y2)), % = 0
% 
%         PercentError(z1,z2), % = 1.3017e-012
%         PercentError(w1,w2), % = 4.3409e-011
% 
% 
% 
% 
% 
% 
%     %%%MEMORY FOOTPRINT
% 
%     whos Matrix Object
% 
% 
%     %       Name           Size                   Bytes  Class       Attributes
%     % 
%     %       Matrix      2400x15000            288000000  double                
%     %       Object      2400x15000                 8102  KronProd    
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USAGE:
%
%  M=KronProd(OPSET,OPINDS,DOMAINSIZES,SCALARCOEFF)
%
%in:
%
% 
%         OPSET  - a cell array containing the operands of the Kronecker product 
%                  as an unordered set.
%  
%                  The order in which the operands will be Kronecker multiplied
%                  is specified by the OPINDS argument (see below). Note that even if some 
%                  of the operands of the Kronecker product are identical to each other, 
%                  only one copy needs to be reside in OPSET. The OPINDS argument
%                  can be used to indicate that the same member of OPSET is to be used
%                  multiple times. Efficiency of various class methods is
%                  increased this way.
%
%                  Note also that a scalar member of OPSET can be used to indicate  
%                  an operand that is multiple of the identity matrix, c*eye(N)  
%                  as opposed to storing a full-sized matrix of this form.
%                  For example, OPSET={A,B,pi} means A,B,and pi*eye(N) will be
%                  Kronecker multiplied with each other in some order. The value of N 
%                  here is determined by the DOMAINSIZES argument (see below).
%
%                  One restriction, however, is that each scalar member of
%                  OPSET can only be used to represent a unique matrix.
%                  Hence, if one desires the Kronecker product of pi*eye(N)
%                  and pi*eye(M), with M~=N, then OPSET={pi,pi} must contain pi twice
%                  and DOMAINSIZES=[M,N].
%
%                  
% 
%         OPINDS - a vector of indices into OPSET. The order of the indices specifies in what
%                  order (from right to left) the matrices/operators in OPSET will appear in the
%                  Kronecker product. For example, OPSET={A,B,C} with OPINDS=[1 2 3]
%                  would correpond to kron(C,kron(B,A)) or equivalently, a tensorial 
%                  transform which multiplies A with all the columns, then multiplies    
%                  B with all the rows, and finally C with all 3rd-dimensional axes of an
%                  appropriately sized 3D array.
%                  Alternatively, OPINDS=[1,1,2] would be equivalent to kron(B,kron(A,A)),
%                  which applies A to both the rows and columns and B to the 3rd dimension.
%                 
%                  Default: the default value of OPINDS is 1:length(OPSET)
%                                             
% 
%         DOMAINSIZES - A vector giving the dimensions of the domain spaces of each operand.                     
%                       The components DOMAINSIZE(i) must be set as follows
%
%                        DOMAINSIZES(i)=size( OPSET{OPINDS(i)} , 2)    (Eq. 2)
%                      
%                       with the following exceptions:
%                        
%                        (1) If OPSET{OPINDS(i)}=c where c is a scalar acting 
%                        as a placeholder for c*eye(N), then one must set
%                        DOMAINSIZES(i)=N.
%                            
%                        (2) For convenience, DOMAINSIZES(i) can be set to
%                        NaN or Inf. This will be treated the same as if
%                        it were set according to Eq. 2 for that particular
%                        index i.
%                        
%                        (3) For convenience, if the entire vector DOMAINSIZES is
%                        set to the empty matrix, to NaN, or to Inf then
%                        this will be treated as if it were set according
%                        to Eq. 2 for all i.
%
%                        The default value of DOMAINSIZES is [].
%
% 
%         SCALARCOEFF - (optional) global scalar coefficient applied to the entire Kronecker
%                       product. Default = +1
%
%
%out:
%
%          M:  The KronProd object. In addition to the various class
%              methods mentioned above, it has the following class properties 
%              all get-accessible by dot indexing syntax:
%
%                 M.opset: The constructor input OPSET described above, but possibly
%                          reordered and with unused operands discarded.
%                 M.opinds: The constructor input OPINDS described above, but possibly
%                           modified to account for reordering of OPSET.
%
%                 M.numops: same as length(opset)
%                 M.eyemask: Boolean vector of length M.numops indicating which members of
%                            M.opset are scalars standing in for matrices c*eye(N).                           
%                 M.domset: Vector of length M.numops indicating
%                           size(M.opset{i},2) with adjustment for members
%                           of M.opset that are scalars standing in for
%                           matrices c*eye(N).
%
%                           
%                 M.scalarcumprods: A deprecated property.
%
%                 M.maxdim: The number of operands in the Kronecker product, length(opinds).                
%                 M.scalarcoeff: The scale factor 's' participating in the product.
%                                It will change appropriately if M is
%                                multiplied by a scalar. 
%                 M.scalarmask:  Boolean vector of length M.maxdim, indicating which 
%                                tensorial dimensions use a scalar to represent a
%                                multiple of an identity matrix, c*eye(N).       
%                 M.domainsizes: A vector such that M.domainsizes(i) is the dimension 
%                                of the domain space of the i-th operand.
%                                Equivalently, if M is to be used as a
%                                tensorial transform, then M.domainsizes is the size
%                                of the array on which it will operate.                                
%                 M.rangesizes = Analogous to M.domainsizes, this vector gives the
%                                the dimensions of the array output from the tensorial
%                                transform applied by M.



properties
  opset={};
  opinds=[];
  numops=0;
  eyemask=[];
  domainset=[];
  maxdim=0;
  scalarcoeff=+1;
  scalarcumprods=1;
  scalarmask=[];

  domainsizes=[];
  rangesizes=[];

end
 

methods
    
    function M=KronProd(opset,opinds,domainsizes,scalarcoeff)

          if ~nargin, return; end

          if ~iscell(opset), opset={opset}; end  

          if nargin<2, opinds=1:length(opset); end
          
          if nargin<3, 
            domainsizes=[];
          end   

         if nargin>=4, 
            M.scalarcoeff=full(double(scalarcoeff));
         end


         %%%Weed out unused operands - stray operands will impact upon efficiency

           [qq,ii,jj]=unique(opinds);
           opset=opset(qq);
           opinds=jj;


         %%%End weed


         M.opset=opset;
         M.opinds=opinds(:).';
         M.numops=length(M.opset);
         M.maxdim=length(M.opinds);
         M.scalarcumprods=ones(size(M.opinds)); 
         M.scalarmask=logical(row0s(M.maxdim));
         %M.domainsizes=domainsizes(:).';
         M.rangesizes=ones(size(M.domainsizes));


         %%process new options for domainset
         %nColSet=cellfun('size',M.opset,2);
         nColSet=cellfun(@(x)size(x,2),M.opset);
         nColMap=nColSet(M.opinds);
         nonfinMap=~isfinite(domainsizes);
         domainsizes=domainsizes(:).';


         if  (numel(domainsizes)==1 && nonfinMap ) || isempty(domainsizes) 
            domainsizes=nColMap; 
         else    
            domainsizes(nonfinMap)=nColMap(nonfinMap);  
         end
         M.domainsizes=domainsizes;

         if length(M.domainsizes)~=M.maxdim,
           error('Inconsistency between domain sizes and operator index set.');
         end


         cc=0;

         for ii=M.opinds

            cc=cc+1;
            sz=size(M.opset{ii});

            if isnum(M.opset{ii})
              M.rangesizes(cc)=M.domainsizes(cc);
              M.scalarmask(cc)=logical(1);
              M.scalarcumprods(cc)=M.opset{ii};
            elseif sz(2)==M.domainsizes(cc)
              M.rangesizes(cc)=sz(1);
            else
             error(['DOMAINSIZE(' num2str(cc) ') not consistent with operand dimensions.']);
            end

         end

         M.scalarcumprods=cumprod(M.scalarcumprods);
         M.eyemask=false(1,M.numops);
         M.domainset= nColSet;     

             if any(M.scalarmask),

                 ScalarMap=[ M.opinds(M.scalarmask) ; M.domainsizes(M.scalarmask) ];
                 ScalarMap=unique(ScalarMap.','rows');
                 opinds_shrink=unique(ScalarMap(:,1));

                 if length(opinds_shrink)<length(ScalarMap(:,1));
                    error 'A scalar member of OPSET is being used to represent two or more differently sized matrices' 
                 else

                      M.eyemask(ScalarMap(:,1))=true;
                      M.domainset(ScalarMap(:,1))=ScalarMap(:,2).';
                 end

             end
 

    end
    
    function M=set.scalarcoeff(M,val)
        
        M.scalarcoeff=full(double(val));
        
    end
 end
end