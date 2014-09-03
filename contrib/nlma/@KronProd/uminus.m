function M=uminus(M)
%Unary minus for KronProd class


M.scalarcoeff=full(double(-M.scalarcoeff));

