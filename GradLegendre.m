function [dP] = GradJacobiP(r,N)

% function [dP] = GradJacobiP(r, alpha, beta, N);
% Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
%	       at points r for order N and returns dP[1:length(r))]        
N=N-1;
dP = zeros(length(r), 1);
if(N == 0)
  dP(:,:) = 0.0; 
else
  dP = sqrt(N*(N+1))*JacobiP(r(:),1,1, N-1);
end
return
