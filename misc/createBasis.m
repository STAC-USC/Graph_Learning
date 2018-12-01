function [U,lambdas] = createBasis(M,order)
[U,eigenval] = eig(M);
lamb = diag(eigenval);
[lambdas,index] = sort(lamb,order);
U = U(:,index);
% change the maximum element of each eigenvector to positive
[~,I] = max(abs(U),[],1);
I = sub2ind(size(U),I,1:size(U,2));
U = bsxfun(@times,U,sign(U(I)));
