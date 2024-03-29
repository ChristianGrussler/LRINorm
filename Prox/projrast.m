%PROJRAST Projection onto the epi-graph of the scaled low-rank inducing 
%   Frobenius and spectral norms for integer-valued r. 
%  
%   [X,xv,final] = PROJRAST(Z,zv,r,p,gamma) computes the projection onto the
%   epi-graph of the scaled low-rank inducing norm ||.||_{\ell_p,r*} 
%   evaluated in (Z,zv) for integer-valued r > 0, i.e. 
%   a) (X,xv) is solution to the optimization problem: 
%       
%       minimize_(X,xv)  ||X-Z||_{\ell_2}^2+(xv-zv)^2
%                   s.t. xv >= gamma*||X||_{\ell_p,r*}
%   where p=2 or p=inf.
%   b) final.t and final.s are the final values of parameters to the two 
%   nested search. For p = inf, final.k is the final value of the third inner 
%   search parameter. 
%   If ||Z/gamma||_{\ell_p^D,r} <= 1, then t = s = k = []. 
%
%   [...] = PROJRAST(Z,zv,r,p,gamma,option) allows us to specify furhter options:
%       1. [...] = PROJRAST(Z,zv,r,p,'vec') is used to flag that the
%       vector-valued problem is to be solved.
%       2. [...] = PROJRAST(Z,zv,r,p,gamma,'tol',tol) sets the relative
%       tolerance of the deciding about zeros e.g. if Z is matrix then
%       for all i: \sigma(X)_i = 0 if |\sigma(X)_i| <= tol ||Z||_{\ell_p,r). 
%       Default value: tol = 1e-12.
%       3. [...] = PROJRAST(Z,zv,r,p,gamma,...,'search',search) changes 
%       from default binary search to linear search over
%       a) t if search.t = 0.
%       b) s if search.s = 0.
%       c) k if search.k = 0 (only for p =inf).
%       4. [...] = PROJRAST(Z,zv,r,p,gamma,...,'init',init) changes 
%       from default binary search start values t_0 = 1, k_0 = 1, s_0 = 0
%       to 
%       a) t_0 if init.t = t_0.
%       b) s_0 if init.s = s_0.
%       c) k_0 if init.k = k_0 (only for p =inf).
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
function [X,xv,final] = projrast(Z,zv,r,p,gamma,varargin)
% Check if p = 2 or p = inf
if p ~= 2 && p~=inf
    error('p can only be equal to 2 or inf');
end

%% Check options
dim = size(Z);
min_mn = min(dim);
max_mn = max(dim);
vec = 0; % Flag Z has matrix
for i = 1:length(varargin)
     
    if strcmp(varargin{i},'vec') 
        if min_mn==1 
            if r <= max_mn
                vec = 1; %Flag Z as vector
            else
                error('r is larger than length(Z)');    
            end
        else
            error('Z is not a vector'); 
        end
    end
   
end

if p == inf
    p = 1;
end

if vec == 1 
    %% Vector valued problem
    % Moreau decomposition
    [Y,w,final] = projrnorm(Z,zv,r,p,gamma,varargin);
    X = Z - Y;
    xv = zv - w;
else
    %% Matrix valued problem
    [U,S,V] = svd(Z);
    d = diag(S);
    [y,w,final] = projrnorm(d,zv,r,p,gamma,varargin);
    % Moreau decomposition
    S(1:min_mn,1:min_mn) = diag(d - y);
    X = U*S*V';
    xv = zv-w;
end

end
