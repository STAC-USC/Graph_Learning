clear 
clc
close all


my_eps_outer = 1e-4; my_eps_inner = 1e-6; max_cycles = 40;
scale = 1;
isNormalized = 0;

load('animals.mat');
% % compute correlation matrix
S = cov(data',1); 
% for binary data we add +1/3 to diagonals(suggested by Banerjee et al. ''Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data (2008)
S = S  + (1/3)*eye(size(S));
A_mask=ones(size(S)) - eye(size(S));
alpha = 0.00;
[Laplacian,~,convergence] = estimate_ggl(S,A_mask,alpha,my_eps_outer,my_eps_inner,max_cycles,2);
%[Laplacian,~,convergence] = estimate_ddgl(S,A_mask,alpha,my_eps_outer,my_eps_inner,max_cycles,2);
%[Laplacian,~,convergence] = estimate_cgl(S,A_mask,alpha,my_eps_outer,my_eps_inner,max_cycles,2);

Laplacian(abs(Laplacian) < my_eps_outer) = 0;  % threshold 
draw_animal_graph(Laplacian,names);