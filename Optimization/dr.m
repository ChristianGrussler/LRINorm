%DR  Douglas-Rachford proximal splitting algorithm with two proximal 
%    mappings as input
%   
%   X = DR(@(Z,gamma)prox_f(Z,gamma,param_f),@(Z,gamma)prox_g(Z,gamma,param_g),dim) 
%   runs the Douglas-Rachford splitting algorithm for the objective function 
%   f + g, where:
%       1. prox_f(z,gamma,param_f) and prox_g(z,gamma,param_g) are funcionts 
%       that derive the proximal mappings of gamma*f and gamma*g, depending
%       on possible parameters param_f and param_g.
%       2. dim is the dimension of the arguments of f and g. 
%   
%   Note: If f and g are convex, then X = argmin(f(X) + g(X)).
%   
%   [X,Z_fix,iter,D] = DR(@(Z,gamma)prox_f(Z,gamma,param_f),@(Z,gamma)prox_g(Z,gamma,param_g),dim) also returns:
%       1. Z_fix = fix point of that Douglas-Rachford iterations
%       2. iter = total number of Douglas-Rachford iterations
%       3. D = Z_fix - X, which is the solution of the dual problem if f
%       and g are convex. 
%
%   [...] = DR(prox_f,prox_g,dim,option) allows to specify further options:
%       1. [...] = DR(@(Z,gamma)prox_f(Z,gamma,param_f),@(Z,gamma)prox_g(Z,gamma,param_g),dim,...,'rho',rho,...) set the step
%       length update of the fix-point update, i.e.,
%       Z_{k+1} = Z_k + rho*(Y_k - X_k), where 0 < rho < 2.
%       The default value is rho = 1.
%       2. [...] = DR(@(Z,gamma)prox_f(Z,gamma,param_f),@(Z,gamma)prox_g(Z,gamma,param_g),dim,...,'gamma',gamma,...) sets gamma 
%       to another value than the default gamma = 1.
%       3. [...] = DR(@(Z,gamma)prox_f(Z,gamma,param_f),@(Z,gamma)prox_g(Z,gamma,param_g),dim,...,'Z0',Z0,...) sets the 
%       initial value of the fix-point iteration.
%       The default choice is Z0 = 0.
%       4. [...] = DR(@(Z,gamma)prox_f(Z,gamma,param_f),@(Z,gamma)prox_g(Z,gamma,param_g),dim,...,'tol',tol,...) 
%       sets the tolerance for zero entries as well as to stop the 
%       iteratrions once norm(Y_k-X_k,'fro') < tol
%       The default tol-value is sqrt(eps).

function [X,Z_fix,iter,D] = dr(prox_f,prox_g,dim,varargin)

% Set default values
rho = 1; % Default value for rho
Z = zeros(dim); % Default value for Z0
tol = sqrt(eps); % Default tolerance
gamma = 1; % Default value for gamma

% Initialize iteration counter
iter = 0;

% Set optional values
for i = 2:2:length(varargin)
   if strcmp(varargin{i-1},'rho')
       rho = varargin{i};
   elseif strcmp(varargin{i-1},'Z0')
       Z = varargin{i};
   elseif strcmp(varargin{i-1},'tol')
       tol = varargin{i};
   elseif strcmp(varargin{i-1},'gamma')
       gamma = varargin{i};
   end
end

  
   
% Initialize X, Y and error between X and Y   
X = zeros(dim);
Y = zeros(dim);
err_XY = tol+1; % Larger than the tol so that loop starts

% Set step-size for printing err_XY
err_bound = 1e-1; % Error bound step-size for display 

% Compute the Douglas-Rachford iterations

while err_XY >= tol

    iter = iter+1; % Increase counter
                 
    % Display error between X_iter and Y_iter
    if err_XY <= err_bound
        disp(['Error between X_', num2str(iter) , ' and ' 'Y_', num2str(iter) ,' <= ', num2str(err_bound)]);
        err_bound = err_bound/10;
    end           
    % Compute Douglas-Rachford steps        
    X = prox_f(Z,gamma);
    Y = prox_g(2*X-Z,gamma);
    Z = Z+rho*(Y-X);
            
    % Update iteration error norm(X_k - Y_k,'fro')
    err_XY = norm(X-Y,'fro');           
end

% Set the final solution
X = Y;

% Compute dual variable D
D = (Z-X)/gamma;
D(abs(D) < tol) = 0;
% Set fix point
if iter > 0
     Z_fix = Z;
else
    % Set fix point if iter = 0 and D = 0
    Z_fix = X; 
end

X(abs(X) < tol) = 0;

end