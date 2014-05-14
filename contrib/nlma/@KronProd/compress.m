function M=compress(M)
%COMPRESS method for KronProd.
%
%Weeds redundant members from M.OPSET

[opset,opinds]=trimobjset(M.opset,M.opinds);

M=KronProd(opset,opinds,M.domainsizes,M.scalarcoeff);
