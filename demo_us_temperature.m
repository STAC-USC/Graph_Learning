clear all
close all

% load us temperature data
load('us_temp_data.mat');

% sample covariance
S = cov(data_2d_state,1);

% initializaions
filter_type = 'exp';
beta = 0.5; % filter parameter
graph_filter_ideal = @(x)(graph_filter_fwd(x,beta,filter_type) );
[U,sigma_sq_C] = createBasis(S,'descend');
max_sigma=(max(sigma_sq_C));
sigma_orig = sigma_sq_C/max_sigma;

% step I: prefilter
sigma_sq_C = sigma_sq_C/max_sigma; sigma_sq_C(sigma_sq_C <= 10^-10) = 0;
lambdas_current = graph_filter_inv(sigma_sq_C,beta,filter_type);
orig_sigmas = 1./lambdas_current; orig_sigmas(orig_sigmas==Inf)=0;
S_prefiltered = U * diag(orig_sigmas) * U';

% step II: graph learning 
Laplacian = estimate_cgl(S_prefiltered,ones(size(S_prefiltered)),0.000,10^-5,10^-7,40);

% step III: filter parameter estimation (for a desired filter type a filter parameter selection step)
%  Note: for exponential filter filter parameter selection step can be skipped, 
%        becayse the output graphs are scaled versions of eachother for different beta parameter
%        please refer to the paper for further details

% show resulting graph on the US map
draw_us_temp_graph(Laplacian, center_vector);
