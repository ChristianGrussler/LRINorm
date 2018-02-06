%PROXNORMRAST Prox of low-rank inducing Frobenius and spectral norm for 
%   integer-valued r
%  
%   [X,final] = PROXNORMRAST(Z,r,p) produces the prox of half of the low-rank 
%   inducing norm ||.||_{\ell_p,r*} evaluated in Z for real-valued r > 0, 
%   i.e. 
%   a) X is solution to the optimization problem: 
%       
%       minimize_X ||X||_{\ell_p,r*} + 0.5*||X-Z||_{\ell_2}^2,
%   where p=1 or p=2
%   b) final.t and final.s are the final values of parameters to the two 
%   nested search. For p = 1, final.k is the final value of the third inner 
%   search parameter. 
%   If ||Z/gamma||_{\ell_p^D,r} <= 1, then t = s = k = []. 
%
%   [...] = PROXNORMRAST(Z,r,p,option) allows to specify furhter options:
%       1. [...] = PROXNORMRAST(z,r,p,'vec') is used to flag that the
%       vector-valued problem is to be solved.
%       2. [...] = PROXNORMRAST(z,r,p,'gamma',gamma) produces the prox of 
%       half of the low-rank inducing norm ||.||_{\ell_p,r*} multiplied by 
%       gamma and evaluated at Z, i.e. X is the solution to the 
%       optimization problem: 
%           minimize_X gamma*||X||_{\ell_p,r*} + 0.5*||X-Z||_{\ell_2}^2,
%       where p=1 or p=2.
%       3. [...] = PROXNORMRAST(z,r,p,'tol',tol) sets the relative
%       tolerance of the deciding about zeros e.g. if Z is matrix then
%       for all i: \sigma(X)_i = 0 if |\sigma(X)_i| <= tol ||Z||_{\ell_p,r). 
%       Default value: tol = 1e-12.
%       2. [...] = PROXNORMRAST(z,r,p,...,'search',search) changes 
%       from default binary search to linear search over
%       a) t if search.t = 0.
%       b) s if search.s = 0.
%       c) k if search.k = 0 (only for p =1).
%       3. [...] = PROXNORMRAST(z,r,p,...,'init',init) changes 
%       from default binary search start values t_0 = 1, k_0 = 1, s_0 = 0
%       to 
%       a) t_0 if init.t = t_0.
%       b) s_0 if init.s = s_0.
%       c) k_0 if init.k = k_0 (only for p =1).
function [X,final] = proxnormrast(Z,r,p,varargin)
%% Check options
[n,m] = size(Z);
min_mn = min([n m]);
max_mn = max([n m]);
vec = 0; % Flag Z has matrix
gamma = 1; % Initialize gamma
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
    
    if strcmp(varargin{i},'gamma') 
        gamma = varargin{i+1};
    end
    
end
    
if vec == 1 
    %% Vector valued problem
    % Moreau decomposition
    [sol,final] = projrnormball(Z/gamma,r,p,varargin);
    X = Z - gamma*sol;
else
    %% Matrix valued problem
    [U,S,V] = svd(Z);
    d = diag(S);
    [sol,final] = projrnormball(d/gamma,r,p,varargin);
    % Moreau decomposition
    S(1:min_mn,1:min_mn) = diag(d - gamma*sol);
    X = U*S*V';
end
end