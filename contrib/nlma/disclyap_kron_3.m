function X=disclyap_kron_3(G,V,symmetric)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% disclyap_kron_3.m
%
% Uses the doubling algorithm to solve the discrete Lyapunov equation 
% X=kron_power(G,j)*X*kron_power(G',k)+V
%
% if symmetric =1, then X and V are assumed symmetric
%
% Alexander Meyer-Gohde 15/08/2012
% 
%
%Included in the package nonlinear_MA_mom
%
%THIS VERSION: 1.0.8 April 23, 2013
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if symmetric==1
    P0=0.5*(triu(V)+tril(V)');P0=triu(P0)+triu(P0,1)';
else
    P0=V;
end
A0=G;
matd=1; 
while matd > 1e-12;
    A0_kron=KronProd({A0,A0,A0,1},[1,2,3], [],1);
    %P1=sparse_kron_prod_3(P0',A0',A0',A0');
    %P1=sparse_kron_prod_3(P1',A0,A0,A0);
    P1=P0+A0_kron*P0*A0_kron';  
    A1=A0*A0; 
    matd=max(max(abs(P1-P0))); 
    P0=P1;A0=A1; 
end 

if symmetric==1
    X=0.5*(triu(P0)+tril(P0)');X=triu(X)+triu(X,1)';
else
    X=P0;
end
