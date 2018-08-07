function [M,rankM,err,D,Z_fix,iter] = drhankelapprox(H,r,varargin)
%DRHANKELAPPROX Douglas-Rachford proximal splitting for low-rank 
%   approximation with Hankel constraint through the low-rank inducing 
%   Frobenius norm and non-convex Douglas-Rachford.
%           
%   M = DRHANKELAPPROX(H,r) determines a Frobenius norm 
%   low-rank Hankel approximation through the Frobenius norm low-rank
%   inducing norm and Douglas-Rachford splitting, i.e., M is the solution 
%   to
%   
%   minimize 0.5*||H||_{\ell_2}^2 - trace(M'H) + 0.5*||M||_{\ell_2,r*}^2 
%        s.t. M is Hankel
%    
%   [M,rankM,e,D,Z_fix,iter] = DRHANKELAPPROX(H,r) also returns:
%       1. rankM = rank(M)
%       2. err = norm(H-M,'fro')
%       3. D = solution of the dual problem
%       4. Z_fix = fix point of the Douglas-Rachford iterations
%       5. iter = total number of Douglas-Rachford iterations
%
%   [...] = DRHANKELAPPROX(H,r,option) allows to specify further options:
%       1. [...] = DRHANKELAPPROX(H,r,...,'solver','NDR',...) changes solver 
%       to use the non-convex Douglas-Rachford.
%       2. [...] = DRHANKELAPPROX(H,r,...,'gamma',gamma,...) multiplies the
%       objective functions with gamma when determining the prox of them.
%       The default value is set to gamma = 1.
%       3. [...] = DRHANKELAPPROX(H,r,...,'rho',rho,...) set the step length
%       update of the fix-point update, i.e. 
%       Z_{k+1} = Z_k + rho*(Y_k - X_k), where 0 < rho < 2.
%       The default value is rho = 1.
%       4. [...] = DRHANKELAPPROX(H,r,...,'Z0',Z0,...) sets the initial value of
%       the fix-point iteration,. The default choice is Z0 = zeros(size(H)).
%       5. [...] = DRHANKELAPPROX(H,r,...,'tol',tol,...) sets the relative
%       tolerance of:
%           + The numerical rank: rankM = rank(M/norm(H,'fro'),tol)
%           + Iterations stop: Stop if (Y_k -X_k)/norm(H,'fro') < tol
%           + Zero values: E.g., D(abs(D/norm(H,'fro')) < tol) = 0
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
dim = size(H);
if r > min(dim)
    error('r is larger than min(size(H))');
end

%% Set default values
solv = 0; % Solve low-rank inducing norm problem
tol = sqrt(eps); % Default tolerance value
gamma = 1; % Default gamma value
rho = 1; % Default rho value
Z0 = zeros(dim); % Default Z0

%% Read of optional values
for i = 2:2:length(varargin)
   if strcmp(varargin{i-1},'solver')
       if strcmp(varargin{i},'NDR')
           solv = 1; % Use non-convex Douglas-Rachford
       end
   elseif  strcmp(varargin{i-1},'tol')
       tol = varargin{i};
   elseif strcmp(varargin{i-1},'gamma')
       gamma = varargin{i};
   elseif strcmp(varargin{i-1},'rho')
       rho = varargin{i};
   elseif strcmp(varargin{i-1},'Z0')
       Z0 = varargin{i};
   end
end

% Define absolute tolerance
norm_H = norm(H,'fro'); % Compute norm(H,'fro')
tol_H = norm_H*tol; % Define absolute tolerance


%% Start Douglas-Rachford iterations
% Choose between convex and non-convex Douglas-Rachford
if  solv == 1 
    [M,Z_fix,iter,D] = dr(@(Z,gamma)proxnonconv_square(Z,r,2,gamma),@(Z,gamma)projhankel(Z,gamma,-H),dim,'tol',tol_H,'gamma',gamma,'Z0',Z0,'rho',rho);

else
    [M,Z_fix,iter,D] = dr(@(Z,gamma)proxnormrast_square(Z,r,2,gamma),@(Z,gamma)projhankel(Z,gamma,-H),dim,'tol',tol_H,'gamma',gamma,'Z0',Z0,'rho',rho);
end
          
% Compute error and rank of the approximation
err = norm(M-H,'fro');
rankM = rank(M,tol_H);
end