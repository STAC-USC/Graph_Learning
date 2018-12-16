function A = laplacianToAdjacency(L,epsilon)

sizeMat = size(L); 
offDiag = ones(sizeMat) - eye(sizeMat); offDiag = offDiag == 1;
A = zeros(sizeMat); A(offDiag) = -L(offDiag);
A(A<epsilon) = 0;

end