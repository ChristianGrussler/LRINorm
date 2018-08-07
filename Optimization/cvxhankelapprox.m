%CVXHANKELAPPROX Low-rank approximation with Hankel constraint by the low-rank 
%          inducing Frobenius norm in CVX for real valued r >= 1.
%          
%   M = CVXHANKELAPPROX(H,r) determines a Frobenius norm low-rank Hankel 
%   approximation through the Frobenius norm low-rank
%   inducing norm with CVX, i.e., M is the solution to
%   
%   minimize 0.5*||H||_{\ell_2}^2 - trace(M'H) + 0.5*||M||_{\ell_2,r*}^2 
%        s.t. M is Hankel
%   for real-valued r >= 1.
%   [M,rankM,err,D,Z,iter] = CVXHANKELAPPROX(H,r) also returns:
%       1. rankM = rank(M)
%       2. err = norm(H-M,'fro')
%       3. D = solution of the dual problem
%    [...] = CVXHANKELAPPROX(H,r,'tol',tol) sets the relative
%       tolerance of:
%           + The numerical rank: rankM = rank(M/norm(H,'fro'),tol)
%           + Iterations stop: Stop if (Y_k -X_k)/norm(H,'fro') < tol
%           + Zero values: E.g., D(abs(D/norm(H,'fro')) < tol) = 0
%           + Multiple r-th singular value in D: No multiplicity if 
%                        (\sigma_r(D)-\sigma_{r+1}(D))/norm(H,'fro') > tol.
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
%   arXiv:1612.03186v1.
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

function [M,rankM,err,D] = cvxhankelapprox(H,r,varargin)

dim = size(H);
if r > min(dim)
    error('r is larger than min(size(H))');
end
r_bar = ceil(r);

% Set default tol values
tol = sqrt(eps);

for i = 2:2:length(varargin)
   if  strcmp(varargin{i-1},'tol')
       tol = varargin{i};
   end
end


%% Test if svd_r(H) is already Hankel for integer r
if (r_bar - r) == 0
    [U,S,V] = svd(H);
    norm_H = norm(diag(S)); % Compute norm(H,'fro')
    M = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Compute low-rank approximation or rank r
    Mh = projhankel(M,zeros(dim),0); % Project onto Hankel subspace
    tol_H = norm_H*tol; % Define absolute tolerance
    Mh(abs(Mh)<tol_H) = 0; %Set small values to zero

    if norm(Mh-M,'fro') >= tol_H
        usecvx = 1; % Flag to start CVX
    else
        usecvx = 0; % Flag to not start CVX
    end
end

%% Start CVX
if usecvx ~= 0
    
    % Define dimensions
    n = dim(1);
    m = dim(2);
    
    % Compute optimal dual variable D
    cvx_begin SDP
    cvx_precision best
        variable T(n,n) symmetric
        variable g(1)
        variable D(n,m)
            minimize 1/2*(trace(T) - g*(n-r))
            subject to   
            [T H+D;
               transpose(H+D) eye(m)] >= 0;
             T >= g*eye(n);
             for i = m-1:-1:-n+1
                sum(diag(fliplr(D),i)) == 0;
             end
    cvx_end    
 

        %% Compute M by the dual formulation
        cvx_begin SDP 
        cvx_precision best
            variable W(n,n) symmetric
            variable P(m,m) symmetric
            variable M(n,m) hankel
                minimize 1/2*trace(P)-trace(H'*M)
                subject to
                [eye(n) - W M; transpose(M) P] >= 0;
                trace(W) == n - r;
                W >= 0;
        cvx_end
    
        
else
    D = H;
end

D(abs(D) < tol_H) = 0; % Set small entries to zero
M(abs(M) < tol_H) = 0; % Set small entries to zero

% Compute err and rank of the approximation
err = norm(M-H,'fro');
rankM = rank(M,tol_H);

end