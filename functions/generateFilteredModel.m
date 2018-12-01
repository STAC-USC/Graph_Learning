function [h_L,h_L_SquareRootMtx] = generateFilteredModel(L,graph_filter)

[U,lambdas_L] = createBasis(L,'ascend'); lambdas_L(lambdas_L<= 10^-10)=0;
sigmas_C = graph_filter(lambdas_L); sigmas_C(sigmas_C <= 10^-10) = 0; 

h_L = U * diag(sigmas_C) * U'; h_L = 0.5*(h_L + h_L');
h_L_SquareRootMtx =  U * diag(sigmas_C.^(0.5));
end
