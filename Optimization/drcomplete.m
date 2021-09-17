%DRCOMPLETE Douglas-Rachford proximal splitting for low-rank 
%   completion through the low-rank inducing Frobenius/spectral norm and 
%   non-convex Douglas-Rachford for integer-valued r.
%      
%   M = DRCOMPLETE(N,Index,r,p) determines a low-rank matrix completion
%   solution though the low-rank inducing Frobenius/spectral norm, i.e., M
%   is a solution to 
%
%   minimize ||M||_{\ell_p,r*}
%       s.t.   M(Index) = N(Index)
%   for p = 2 or p = inf, where Index is the logical incident matrix of the 
%   known entries in N. Therefore, N should contain the correct known entries,
%   whereas the unknowns can be arbitrary.
%
%   [M,rankM,err,D,Z_fix,iter] = DRCOMPLETE(N,Index,r,p) also returns:
%       1. rankM = rank(M)
%       2. err = norm(N-M,'fro') (This is only to compare the performance) 
%       with other methods, when N is explicitly known. 
%       3. D = solution of the dual problem
%       4. Z_fix = fix point of the Douglas-Rachford iterations
%       5. iter = total number of Douglas-Rachford iterations
%
%   [...] = DRCOMPLETE(N,Index,r,p,option) allows to specify further options:
%       1. [...] = DRCOMPLETE(N,Index,r,p,...,'solver','NDR',...) sets the 
%       solver to the non-convex Douglas-Rachford.
%       2. [...] = DRCOMPLETE(N,Index,r,p,...,'gamma',gamma,...) multiplies 
%       the objective functions with gamma when determining the prox of them.
%       The default value is set to gamma = 1.
%       3. [...] = DRCOMPLETE(N,Index,r,p,...,'rho',rho,...) set the step 
%       length update of the fix-point update, i.e. 
%       Z_{k+1} = Z_k + rho*(Y_k - X_k). The default value is rho = 1.
%       4. [...] = DRCOMPLETE(N,Index,r,p,...,'Z0',Z0,...) sets the initial 
%       value of the fix-point iteration, i.e. Z_0 = Z0. The default choice 
%       is Z0 = randn(size(N)).
%       5. [...] = DRCOMPLETE(N,Index,r,p,...,'tol',tol,...) sets the relative
%       tolerance of:
%           + The numerical rank: rankM = rank(M/norm(N0,'fro'),tol)
%           + Iterations stop: Stop if (Y_k -X_k)/norm(N0,'fro') < tol
%           + Zero values: E.g., D(abs(D/norm(N0,'fro')) < tol) = 0
%       where N0(Index) = N(Index) and zero otherwise. 
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
function [M,rankM,err,D,Z_fix,iter] = drcomplete(N,Index,r,p,varargin)

dim = size(N);
if r > min(dim)
    error('r is larger than min(size(N))');
end
% Set default values
solv = 0;
gamma = 1;
rho = 1;
Z0 = randn(dim);
tol = sqrt(eps);

% Set optional values
for i = 2:2:length(varargin)
   if strcmp(varargin{i-1},'solver')
       if strcmp(varargin{i},'NDR')
           solv = 1;
       end
   elseif strcmp(varargin{i-1},'gamma')
       gamma = varargin{i};
   elseif strcmp(varargin{i-1},'rho')
       rho = varargin{i};
   elseif strcmp(varargin{i-1},'Z0')
       Z0 = varargin{i};
   elseif  strcmp(varargin{i-1},'tol')
       tol = varargin{i};
   end
end

% Make Index logical
Index = logical(Index);

% Define matrix with known entries of N and zeros otherwise
N0 = zeros(dim);
N0(Index) = N(Index);

% Set absolut tolerance 
tol_N0 = norm(N0,'fro')*tol;

%% Start Douglas-Rachford iterations

% Choose between convex and non-convex Douglas-Rachford
if  solv == 1
    if p == 2
        [M,Z_fix,iter,D] = dr(@(Z,gamma)proxnonconv(Z,r,2,gamma),@(Z,gamma)projindex(Z,N,Index,zeros(dim),gamma),dim,'tol',tol_N0,'gamma',gamma,'Z0',Z0,'rho',rho);
    else % p == inf
        [M,Z_fix,iter,D] = dr(@(Z,gamma)proxnonconv(Z,r,inf,gamma),@(Z,gamma)projindex(Z,N,Index,zeros(dim),gamma),dim,'tol',tol_N0,'gamma',gamma,'Z0',Z0,'rho',rho);
    end
elseif p == 2
    [M,Z_fix,iter,D] = dr(@(Z,gamma)proxnormrast(Z,r,2,gamma),@(Z,gamma)projindex(Z,N,Index,zeros(dim),gamma),dim,'tol',tol_N0,'gamma',gamma,'Z0',Z0,'rho',rho);
else % p == inf
    [M,Z_fix,iter,D] = dr(@(Z,gamma)proxnormrast(Z,r,inf,gamma),@(Z,gamma)projindex(Z,N,Index,zeros(dim),gamma),dim,'tol',tol_N0,'gamma',gamma,'Z0',Z0,'rho',rho);
end

% Compute error and rank of the approximation
err = norm(M-N,'fro');
rankM = rank(M,tol_N0);
end
