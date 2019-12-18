function [P] = JacobiP(x,N)

%function [P] = JacobiP(x,alpha,beta,N);
% Purpose: Calculate the Legendre function using Jacobi function, at alpha=0 and beta=0;

P = JacobiP(x,0,0, N-1);
return