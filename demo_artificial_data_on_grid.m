clear
clc
close all

%% Simulation parameters
my_eps_outer = 10^(-4); my_eps_inner = 10^(-6);
max_iterations = 20;

%% Generate model
% generate random 8x8 grid graph
BS = 8;
load(['Grid_' num2str(BS) 'x' num2str(BS)  '.mat']);
A_connectivity = grid_Con4; A_mask = A_connectivity;
[Ltrue,~,~] = generateRandomGraphFromConnectivity(A_connectivity,3,0,'uniform');
num_vertices = size(A_connectivity,2);

% generate filter
filter_type = 'hop'; beta = 3; % beta-hop filter with true beta parameter is set here
graph_filter_ideal = @(x)(graph_filter_fwd(x,beta,filter_type) );

% generate graph system
[h_L,h_L_sqrMatrix] = generateFilteredModel(Ltrue,graph_filter_ideal);


%% Generate data samples based on h_L
S_data_cell = generateRandomDataFromSqrtM(h_L_sqrMatrix,30); % creates data with 30 samples per vertex
S_data = S_data_cell{1}; S_data = 0.5*(S_data + S_data');


%% Algorithm implementation 
beta_current = 1; % initialize beta
% Eigendecompose data
[U,sigma_sq_C] = createBasis(S_data,'descend');
sigma_sq_C(sigma_sq_C <= 10^-10) = 0;
for repeat=1:max_iterations
 disp(['-- Iteration ' num2str(repeat) '--']);
% Step I: Prefiltering step
lambdas_current = graph_filter_inv(sigma_sq_C,beta_current,filter_type);
current_sigmas = 1./lambdas_current; 
current_sigmas(current_sigmas==Inf)=0;
% construct unfiltered S
S_prefiltered = U * diag(abs(current_sigmas)) * U'; 
S_prefiltered = 0.5*(S_prefiltered + S_prefiltered'); % symmetrize (in case of numerical inaccuracies)
max_eig_of_S_data = max(current_sigmas);
S_prefiltered = S_prefiltered/max_eig_of_S_data; % normalize

% Step II: Graph learning step
Laplacian = estimate_cgl(S_prefiltered,ones(num_vertices),eps,my_eps_outer,my_eps_inner,max_iterations);
Laplacian = Laplacian/max_eig_of_S_data; % normalize

% Step III: Parameter estimation step
estimated_lambdas = diag(U' * Laplacian * U); estimated_lambdas(estimated_lambdas <= 10^-10) = 0;
estimated_sigmas = graph_filter_fwd(estimated_lambdas,repeat,filter_type);
h_L_current = U * diag(estimated_sigmas) * U'; h_L_current = 0.5*(h_L_current + h_L_current');
error_est(repeat) = norm(S_data-h_L_current,'fro')/norm(S_data,'fro');
% convergence criterion
if repeat > 1 && abs(error_est(repeat-1) - error_est(repeat)) > 0.2
     disp('*** Algorithm has converged ***');
     break;
end
beta_current = beta_current + 1;
end


disp(['Estimated beta = ' num2str(beta_current) ]);
disp(' Figure 1 shows the ground truth graph');
disp(' Figure 2 shows the estimated graph');

f1=figure(1);
draw_grid_graph(laplacianToAdjacency(Ltrue,eps),BS);
axis square
movegui(f1,'west');

f2=figure(2);
draw_grid_graph(laplacianToAdjacency(Laplacian,eps),BS);
axis square
movegui(f2,'east');



    