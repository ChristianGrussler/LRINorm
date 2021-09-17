%PROXNONCONV Non-convex prox of the Frobenius and spectral norm + indicator 
%   function of a matrix with at most rank r.
%  
%   X = PROXNONCONV(Z,r,p,gamma) determines the prox of 
%   gamma*||X||_{\ell_p}+i_{rank(X)<=r}(X), where i_{rank(X)<=r} is the
%   idicator function of all matrices with at most rank r. This is, X is a 
%   solution to the optimization problem: 
%       
%       minimize_X gamma*||X||_{\ell_p}+0.5*||X-Z||_{\ell_2}^2
%             s.t. rank(X) <= r
%
%   where p=2 or p=inf.
%  
%   [...] = PROXNONCONV(Z,r,p,gamma,option) allows us to specify furhter options:
%       1. [...] = PROXNONCONV(Z,r,p,gamma,...,'vec') is used to flag that the
%       vector-valued problem is to be solved.
%       2. [...] = PROXNONCONV(Z,r,p,gamma,...,'tol',tol) sets the relative
%       tolerance of the deciding about zeros e.g. if Z is matrix then
%       for all i: \sigma(X)_i = 0 if |\sigma(X)_i| <= tol ||Z||_{\ell_p,r). 
%       Default value: tol = 1e-12.
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
function X = proxnonconv(Z,r,p,gamma,varargin)
% Check if p = 2 or p = inf
if p ~= 2 && p~=inf
    error('p can only be equal to 2 or inf');
end

%% Check options
dim = size(Z);
min_mn = min(dim);
max_mn = max(dim);

for i = 1:length(varargin)
     
    if strcmp(varargin{i},'vec') 
        if min_mn==1 
            if r > max_mn
                error('r is larger than length(Z)');    
            end
        else
            error('Z is not a vector'); 
        end
    end
    
end

% Compute rank-r approximation of Z
[U,S,V] = svd(Z);
Z = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
X = proxnormrast(Z,r,p,gamma,varargin);


end
