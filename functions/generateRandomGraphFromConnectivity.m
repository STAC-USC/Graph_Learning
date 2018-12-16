% Function for generating random data given
%
% (C) Hilmi Egilmez May 2016

function [Ltrue,S_true,S_sqrt_matrix] = generateRandomGraphFromConnectivity(A_connectivity,edge_weight_par,vertex_weight_par,type_weighting)

if nargin < 4
    type_weighting = 'uniform';
    if nargin < 3
        vertex_weight_par = edge_weight_par;
    end
end

num_vertices = size(A_connectivity,2);

random_Seed = 'shuffle';
num_edges = sum(A_connectivity(:))/2;
rng(random_Seed); % seed random number generator
% apply random weights depending on weight_var
if isequal(type_weighting,'gaussian')
    % Gaussian random
    weights = abs(randn(1,num_edges)*sqrt(edge_weight_par))+0.1; %0.01 + exprnd(weight_var,1,num_edges); % abs(randn(1,num_edges)*sqrt(weight_var));
    vertex_weights = abs(randn(1,num_vertices)*sqrt(vertex_weight_par))+0.1; % 0.01 + exprnd(weight_var,1,num_vertices); %abs(randn(1,num_vertices)*sqrt(weight_var));
elseif isequal(type_weighting,'uniform')
    %  Uniform
    weights = (rand(1,num_edges)*(edge_weight_par-0.1))+0.1; %0.01 + exprnd(weight_var,1,num_edges); % abs(randn(1,num_edges)*sqrt(weight_var));
    vertex_weights = (rand(1,num_vertices)*(vertex_weight_par-0.1))+0.1; % 0.01 + exprnd(weight_var,1,num_vertices); %abs(randn(1,num_vertices)*sqrt(weight_var));
else
    error('Wrong input: type_weighting')
end
% create Laplacian
B = convertAdjToIncidence(A_connectivity);
if vertex_weight_par == 0
    Ltrue = B * diag(weights) * B';
else
    Ltrue = B * diag(weights) * B' + diag(vertex_weights);
end
% Create DataSamples
[U,E] = eig(Ltrue);

% scale graph values when it has too large/small e.vals
d_e = diag(E); 

%%% added to avoid determinant < 1 case for vassilis algorithm
d_e(d_e < 10^-10)=1;
det_L = prod(d_e);
if det_L < 1
    E = E*(det_L^(-1/num_vertices)+eps);
    d_e = d_e*(det_L^(-1/num_vertices)+eps);
end

emin = min(d_e); emax = max(d_e);
scaling = 1;
if det_L == Inf || det_L==0
    scales = linspace(emin,emax,10);
    for i=1:length(scales)
        s = scales(i);
        det_L = prod(d_e./s);
        if det_L ~= Inf && det_L~=0
            scaling = s;
            break;
        end
    end
end

if det_L == Inf || det_L == 0
    error(['Function generateRandomGraphFromConnectivity :  Determinant of Laplacian is equal to ' ...
        num2str(det_L) ' -- try different parameters to generate random graph weights' ]);
end

Ltrue = Ltrue/scaling;
E = E/scaling;
E(E<=10^(-10)) = 0;
idx = E>10^(-10);
E(idx) = abs(1./E(idx));
S_true = U* E *U'; 
%%% outputs
S_sqrt_matrix = U*(E.^(0.5));
S_true = 0.5*S_true + 0.5*S_true'; % symmetrize 
Ltrue = 0.5*Ltrue + 0.5*Ltrue';% symmetrize 

end