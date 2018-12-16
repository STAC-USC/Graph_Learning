% Function for generating random data given square root matrix assoc. with graph
%
% (C) Hilmi Egilmez

function [S_cell,D_cell] = generateRandomDataFromSqrtM(SquareRootMatrix,k_over_n_list)
num_covars = length(k_over_n_list);
S_cell  = cell(1,num_covars); % cell of covariance matrices
D_cell  = cell(1,num_covars); % cell of sample distance matrices
num_vertices = size(SquareRootMatrix,2);
degree_of_freedom = num_vertices; %num_vertices*(num_vertices-1)/2;
% degree_of_freedom = num_vertices*(num_vertices-1)/2;
nSamplesList =  int32(round(degree_of_freedom*k_over_n_list)); % number of data samples

y_samples = SquareRootMatrix*randn(num_vertices,max(nSamplesList));

for i=1:length(k_over_n_list)
    y_selected = y_samples(:,1:nSamplesList(i));
    S_temp = cov(y_selected',1);
    S_cell{i} = S_temp;
end

end

