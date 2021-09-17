%CVXCOMPLETE Matirx completion with low-rank inducing Frobenius/spectral
%   norm and CVX.
%
%   M = CVXCOMPLETE(N,Index,r,p) determines the solution of the
%   convexifications of the low-rank inducing Frobenius/spectral norm
%   with CVX,i.e., M is the solution to
%
%       minimize  1/2*||M||_{\ell_p,r*}^2
%           s.t.  M(Index) = N(Index)
%
%   for p = 2 and p=inf, where Index is the logical incident matrix of the 
%   known entries in N. Therefore, N should contain the correct known entries,
%   whereas the unknowns can be arbitrary.
%
%   [M,rankM,err,D,Z,iter] = CVXCOMPLETE(N,Index,r) also returns:
%       1. rankM = rank(M)
%       2. err = norm(N-M,'fro'); This is only to compare the performance
%       with other methods, when N is explicitly known.
%       3. D = solution of the dual problem
%   [...] = CVXCOMPLETE(N,Index,r,'tol',tol) sets the relative
%       tolerance of:
%           + The numerical rank: rankM = rank(M/norm(N(Index),'fro'),tol)
%           + Iterations stop: Stop if (Y_k -X_k)/norm(N(Index),'fro') < tol
%           + Zero values: E.g., D(abs(D/norm(N(Index),'fro')) < tol) = 0
%           + Multiple r-th singular value in D: No multiplicity if
%                        (\sigma_r(D)-\sigma_{r+1}(D))/norm(N(Index),'fro') > tol.
%       The default tol-value is sqrt(eps).
%
%%%%%%%%%%%%%
% References:
%   - C. Grussler and A. Rantzer and P. Giselsson (2018): 
%   "Low-Rank Optimization with Convex Constraints", 
%   IEEE Transactions on Automatic Control, DOI: 10.1109/TAC.2018.2813009.
%
%   - C. Grussler and P. Giselsson (2016):
%   "Low-Rank Inducing Norms With Optimality Interpreations", 
%   SIAM J. Optim., 28(4), pp. 3057–3078.
%
%   - C. Grussler and P. Giselsson (2017):
%   "Local convergence of proximal splitting methods for rank constrained
%   problems", pp. 702-708, IEEE 56th Annual Conference on Decision and Control
%   (CDC), DOI: 10.1109/CDC.2017.8263743.
%
%   - C. Grussler (2017):
%   "Rank reduction with convex constraints", PhD Thesis, 
%   Department of Automatic Control, Lund Institute of Technology, 
%   Lund University, ISBN 978-91-7753-081-7.
%%%%%%%%%%%%%

function [M,rankM,err,D] = cvxcomplete(N,Index,r,p,varargin)

% Check if p = 2 or p = inf
if p ~= 2 && p~=inf
    error('p can only be equal to 2 or inf');
end

dim = size(N);
if r > min(dim)
    error('r is larger than min(size(N))');
end

% Make Index logical
Index = logical(Index);

% Define matrix with known entries of N and zeros otherwise
N0 = zeros(dim);
N0(Index) = N(Index);

% Set default tol values
tol = sqrt(eps);


for i = 2:2:length(varargin)
    if  strcmp(varargin{i-1},'tol')
        tol = varargin{i};
    end
end


% Define absolute tolerance 
tol_N0 = norm(N0,'fro')*tol;


% Reduce computational cost if m < n
if dim(2) < dim(1)
    transp = 1;
    N = N';
    dim = sort(dim,'descend');
else
    transp = 0;
end
n = dim(1);
m = dim(2);

%% Compute optimal dual variable D

% Define index set of the unknown entries in N
Index_not = logical(ones(dim)-double(Index));
if p == 2
    cvx_begin SDP
    cvx_precision best
    variable T(n,n) symmetric
    variable g(1)
    variable D(n,m)
    minimize 1/2*(trace(T) - g*(n-r))-trace(N0'*D)
    subject to
    [T D;
        transpose(D) eye(m)] >= 0
    T >= g*eye(n)
    D(Index_not) == 0;
    cvx_end
    D(abs(D) < tol_N0) = 0; % Set small entries to zero
    
    %% Compute M by the dual formulation
    cvx_begin SDP
    cvx_precision best
    variable W(n,n) symmetric
    variable P(m,m) symmetric
    variable M(n,m)
    minimize 1/2*trace(P)
    subject to
    [eye(n) - W M; transpose(M) P] >= 0;
    W >= 0;
    trace(W) == n - r;
    M(Index) == N0(Index);
    cvx_end
    
else % p == inf
    cvx_begin SDP
    cvx_precision best
    variable T1(n,n) symmetric
    variable T2(m,m) symmetric
    variable g(1)
    variable D(n,m)
    minimize 1/2*square_pos(1/2*(trace(T1)+trace(T2)-(n+m-2*r)*g))-trace(N0'*D)
    subject to
    [T1 D;
        transpose(D) T2] >= 0
    T1 >= g*eye(n)
    T2 >= g*eye(m)
    D(Index_not) == 0;
    cvx_end
    
    cvx_precision best
    cvx_begin SDP
    variable W1(n,n) symmetric
    variable W2(m,m) symmetric
    variable k(1)
    variable M(n,m)
    minimize 0.5*k^2
    
    [k*eye(n)-W1 M;
        transpose(M) k*eye(m)-W2] >= 0
    trace(W1)+trace(W2) == k*(n+m-2*r)
    W1 >= 0
    W2 >= 0
    M(Index) == N(Index)
    cvx_end
    
end

% Set small values to zero
M(abs(M) < tol_N0) = 0;

% Compute error and rank of the approximation
err = norm(M-N,'fro');
rankM = rank(M,tol_N0);

% Transpose back if necessary
if transp == 1
    M = M';
    D = D';
end

end
