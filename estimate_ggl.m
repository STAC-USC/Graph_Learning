% Function for finding the best generalized Laplacian estimate under connectivity constraints
% This function returns the optimal solution of the GGL problem in [1].
%    Reference:
%    [1] H. E. Egilmez, E. Pavez, and A. Ortega, ?Graph learning from data under structural and Laplacian constraints,?
%       CoRR, vol. abs/1611.05181v1,2016. [Online]. Available: https://arxiv.org/abs/1611.05181v1
%
% Inputs: S - sample covariance/input matrix
%         A_mask - 0/1 adjacency matrix where 1s represent connectivity
%         alpha -  regularization parameter
%         prob_tol - error tolerance for the optimization problem (stopping criterion)
%         inner_tol - error torerance for the nonnegative QP solver
%         max_cycle - max. number of cycles
%         regularization_type - two options:
%                               1 (default): | O |_1 (l_1 regularization)
%                               2          : | O |_1,off (l_1 regularization on off-diagonal entries)
%
% Outputs:
% O   -   estimated generalized graph Laplacian (Theta)
% C   -   estimated inverse of O
% convergence  - is the structure with fields 'frobNorm', 'time':
%                'frobNorm': vector of normalized differences between esimtates of O between cycles: |O(t-1)-O(t)|_F/|O(t-1)|_F
%                'time':  runtime of each cycle
%
%
%  (C) Hilmi Enes Egilmez

function [O, C, convergence] = estimate_ggl(S, A_mask, alpha, prob_tol,inner_tol,max_cycle,regularization_type)

%%% set default parameters
if  nargin < 7
    regularization_type = 1;
    if  nargin < 6
        max_cycle = 20;
        if nargin < 5
            inner_tol = 1E-6;
            if nargin < 4
                prob_tol = 1E-4;
                if nargin < 3
                    alpha = 0; % regularization
                end
            end
        end
    end
end
%%% variables
n = size(S, 2);
if regularization_type == 1
    H_alpha = (alpha)*(2*eye(n) - ones(n));
elseif regularization_type == 2
    H_alpha = (alpha)*(eye(n) - ones(n));
else
    error('"regularization_type" can be either 1 or 2.');
end
K = S + H_alpha;

%%% starting value;
O_init = diag((1 ./ diag(K)));     % S \ eye(p);
C = diag(diag(K));
O = O_init;

% Best output
O_best = O; C_best = C;


%%%
frob_norm = [];
time_counter = [];
converged = false;
cycle = 0;
%

while ~converged && cycle < max_cycle
    
    O_old = O;
    tic;
    %%% Inner loop
    for u=1:n
        
        minus_u = setdiff(1:n,u); % index of u complement
        
        % input matrix variables
        k_uu = K(u,u);
        k_u = K(minus_u,u);
        
        % update Ou_inv
        c_u = C(minus_u, u);
        c_uu = C(u,u);
        Ou_i = C(minus_u, minus_u) - (c_u * c_u'./ c_uu);
        
        %%% block-descent variables
        beta = zeros(n-1,1);
        ind_nz = A_mask(minus_u,u) == 1; % non-zero indices
        A_nnls = Ou_i(ind_nz,ind_nz);
        b_nnls = k_u(ind_nz)/k_uu;
        %%% block-descent step
        out = nonnegative_qp_solver(A_nnls, b_nnls,inner_tol);
        beta_nnls = -out.xopt; %%% sign flip
        beta(ind_nz) = beta_nnls;
        o_u = beta;
        o_uu = (1/k_uu) + o_u' * Ou_i * o_u;
        %%% Update the current Theta
        O(u,u) = o_uu;
        O(minus_u, u) = o_u;
        O(u, minus_u) = o_u;
        %%% Update the current Theta inverse
        cu = (Ou_i * o_u)*(k_uu);
        cuu = 1./(o_uu - (o_u' * Ou_i *o_u));
        C(u,u) = cuu; %C(u,u) = k_uu;
        C(u, minus_u) = -cu;
        C(minus_u, u) = -cu;
        %%% use Sherman-Woodbury
        C(minus_u, minus_u) = (Ou_i + ((cu * cu')./ (k_uu)));
        
    end
    time_interval=toc;
    O_best = O; C_best = C;
    cycle = cycle + 1;
    
    %%% time
    time_counter(cycle) = time_interval;
    %%% calculate frob norms
    frob_norm = [frob_norm norm(O_old - O,'fro')/norm(O_old,'fro')];
    
    fprintf('Cycle: %2d      |Theta(t-1)-Theta(t)|_F/|Theta(t-1)|_F = %5.8f \n',cycle, (frob_norm(end)) );
    
    if cycle > 5
        % convergence criterions (based on theta)
        if(((frob_norm(end))) < prob_tol)
            converged = true;
            O_best = O; C_best= C;
        end
    end
end

O = O_best; C = C_best;

convergence.frobNorm = frob_norm;
convergence.time = time_counter;


