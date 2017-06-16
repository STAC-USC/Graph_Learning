
% Nonnegative quadratic program solver based on block pivoting method (Portugal et al., 1994)

% References:
%  [1] L. F. Portugal, J. J. Judice, and L. N. Vicente, ?A comparison of block pivoting and interior-point algorithms 
%      for linear least squares problems  with nonnegative variables,? Math. Comput., vol. 63, no. 208, pp. 625?643, Oct. 1994.
%  [2] M. Slawski, "Problem-specific analysis of non-negative least squares  solvers with a focus on instances with sparse solutions", 2013 (preprint).
%

% This code is a modified version of the solver originally implemented as
% part of nnlslab MATLAB package available online: https://sites.google.com/site/slawskimartin/nnlslab.zip?attredirects=0&d=1

function[out] = nonnegative_qp_solver(A, b, inner_tol)


% initialization
[n, p] = size(A);
AtA = A;
Atb = b;
qbar = 10;
mv = 1;

maxiter = 50*n;
F = false(p, 1);
Lambda = true(p, 1);
lambda = -Atb;
x = zeros(p, 1);
q = qbar;
ninf = p + 1;
iter = 0;
global f;
f = @(x) x' * (AtA * x - 2 * Atb);

flag = 0;

check = check_kkt_cond(lambda, x);

tol = inner_tol;

if check < tol
    out.xopt = x;
    out.check = check;
end

while (check > tol) && (iter < maxiter)
    fH1 = x < -eps;
    H1 = fH1 & F;
    
    fH2 = lambda < -eps;
    H2 = fH2 & Lambda;
    H1H2 = H1 | H2;
    cardH1H2 = sum(H1H2);
    if cardH1H2 < ninf
        ninf = cardH1H2;
        q = qbar;
    else
        if q > 0
            q = q - 1;
        else
            if sum(F) >= min(n,p)
                r = find(H1, 1, 'last');
            else
                r = find(H1H2, 1, 'last');
            end
            if H1(r)
                H1 = false(p, 1);
                H1(r) = true;
                H2  = false(p, 1);
                fH1 = false(p, 1);
                fH1(r) = true; %%find(F == r);
                fH2 = false(p, 1);
            else
                H1  = false(p, 1);
                H2  = false(p, 1);
                H2(r) = true;
                fH1 = false(p, 1);
                fH2 = false(p, 1);
                fH2(r) = true;
            end
        end
    end
    F(fH1) = false;
    F(H2) = true;
    
    Lambda(H1) = true;
    Lambda(H2) = false;
    lfH2 = length(fH2);
    
    lfH1 = length(fH1);
    
    %%% 0.1: rule of thumb.
    if  flag && (lfH1 + lfH2) <= max(1, 0.1 * length(F))
        % use up- and downdating
        % downdating, maintaining ordering
        if lfH1 > 0
            for j=1:lfH1
                dd = find(fH1, j, 'last');
                R = choldownmatlab(R, dd(j));
            end
        end
        % then updating
        if lfH2 > 0
            for j=1:lfH2
                R = cholinsertgram(R, AtA(H2(j), H2(j)), AtA(H2(j), Fprime));
                Fprime(H2(j)) = 1;
            end
        end
    else
        R = chol(AtA(F,F));
        flag = 1;
    end
    x = zeros(p, 1);
    x(F) =  R \ (R' \ Atb(F));
    lambda(Lambda) = AtA(Lambda, F) * x(F) - reshape(Atb(Lambda), sum(Lambda), 1);
    lambda(F) = 0.9 * eps;
    iter = iter + 1;
    check = check_kkt_cond(lambda, x);
end

out.xopt = x;
out.check = check;

end



function[kktopt] = check_kkt_cond(grad, x)

num = abs(grad) .* (grad > eps);
score = num;

Pc = num < eps;

if sum(Pc) > 0
    score(Pc) = num(Pc);
end

if sum(~Pc) > 0
    score(~Pc) = num(~Pc) ./ x(~Pc);
end

Pc = num < eps | score < 1;
P  = ~Pc;

if sum(P) > 0
    kktP = max(abs(x(P)));
else
    kktP = 0;
end

if sum(Pc) > 0
    kktPc = max(abs(grad(Pc)));
else
    kktPc = 0;
end

N = x < -eps;

if sum(N) > 0
    kktN = max(abs(x(N)));
else
    kktN = 0;
end

kktopt = max([kktP kktPc kktN]);

end

% Fast Cholesky insert and remove functions
function R = cholinsertgram(R, diag_k, col_k)
if isempty(R)
    R = sqrt(diag_k);
else
    %  col_k = x'*X; %
    R_k = R'\col_k'; % R'R_k = (X'X)_k, solve for R_k
    R_kk = sqrt(diag_k - R_k'*R_k); % norm(x'x) = norm(R'*R), find last element by exclusion
    R = [R R_k; [zeros(1,size(R,2)) R_kk]]; % update R
end
end



function Lkt =  choldownmatlab(Lt, k)
% This is for non-sparse matrix
% cholesky downdating
% A in R^(n,p)
% G = A'* A = L * L', where L, L' come from cholesky decomposition
% now  removes kth column from A, denoted by Ak. Gk := Ak' * Ak
% Given L' and k, choldown computes the chol. decomposition of  Gk
% i.e. Lk' * Lk = Gk, without processing of A, G

p = length(Lt);

% drop the kth clm of Lt
Temp = Lt;
Temp(:,k) = []; % Temp in R^(p,p-1)

% Givens Rotations
for i = k:p-1,
    a = Temp(i,i);
    b = Temp(i+1,i);
    r = sqrt(sum(Lt(:,i+1).^2) - sum(Temp(1:i-1,i).^2));
    c =  r * a / (a^2+b^2);
    s =  r * b / (a^2+b^2);
    % ith row of rotation matrix H
    Hrowi = zeros(1,p); Hrowi(i) = c; Hrowi(i+1) = s;
    % (i+1)th row of ration matrix H
    Hrowi1 = zeros(1,p); Hrowi1(i) = -s; Hrowi1(i+1) = c;
    % modify the ith and (i+1)th rows of Temp
    v = zeros(2,p-1);
    v(1,i:p-1) = Hrowi * Temp(:,i:p-1);
    v(2,i+1:p-1) = Hrowi1 * Temp(:,i+1:p-1);
    Temp(i:i+1,:) =  v;
end

% drop the last row
Lkt = Temp(1:p-1,:);
end