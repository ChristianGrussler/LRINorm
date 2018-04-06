%PROXNORMRAST_SQUARE Prox of half of the squared low-rank inducing Frobenius 
%   and spectral norms for integer-valued r.
%  
%   [X,final] = PROXNORMRAST_SQUARE(Z,r,p,gamma) produces the prox of half of the low-rank 
%   inducing norm ||.||_{\ell_p,r*} evaluated in Z for integer-valued r > 0, 
%   i.e., 
%   a) X is solution to the optimization problem: 
%       
%       minimize_X 0.5*gamma*||X||_{\ell_p,r*}^2 + 0.5*||X-Z||_{\ell_2}^2,
%
%   where p=2 or p=inf.
%   b) final.t and final.s are the final values of parameters to the two 
%   nested search. For p = inf, final.k is the final value of the third inner 
%   search parameter. 
%  
%   [...] = PROXNORMRAST_SQUARE(Z,r,p,gamma,option) allows us to specify furhter options:
%       1. [...] = PROXNORMRAST_SQUARE(Z,r,p,gamma,...,'vec') is used to flag that the
%       vector-valued problem is to be solved.
%       2. [...] = PROXNORMRAST_SQUARE(Z,r,p,gamma,...,'tol',tol) sets the relative
%       tolerance of the deciding about zeros e.g. if Z is matrix then
%       for all i: \sigma(X)_i = 0 if |\sigma(X)_i| <= tol ||Z||_{\ell_p,r). 
%       Default value: tol = 1e-12.
%       3. [...] = PROXNORMRAST_SQUARE(Z,r,p,gamma,...,'search',search) changes 
%       from default binary search to linear search over
%       a) t if search.t = 0.
%       b) s if search.s = 0.
%       c) k if search.k = 0 (only for p =inf).
%       4. [...] = PROXNORMRAST_SQUARE(Z,r,p,gamma,...,'init',init) changes 
%       from default binary search start values t_0 = 1, k_0 = 1, s_0 = 0
%       to 
%       a) t_0 if init.t = t_0.
%       b) s_0 if init.s = s_0.
%       c) k_0 if init.k = k_0 (only for p =inf).
function [X,final] = proxnormrast_square(Z,r,p,gamma,varargin)
% Check if p = 2 or p = inf
if p ~= 2 && p~=inf
    error('p can only be equal to 2 or inf');
end

%% Check options
[n,m] = size(Z);
min_mn = min([n m]);
max_mn = max([n m]);
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

      
if vec == 1 
    %% Vector valued problem
    % Moreau decomposition
    [Y,final] = projrast(Z,0,r,p,sqrt(gamma),varargin);
    X = Y;
else
    %% Matrix valued problem
    [U,S,V] = svd(Z);
    d = diag(S);
    [y,final] = projrast(d,0,r,p,sqrt(gamma),varargin);
    % Moreau decomposition
    S(1:min_mn,1:min_mn) = diag(y);
    X = U*S*V';
end

end