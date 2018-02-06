% PROJRNORMBALL Projection onto the unit-ball of the truncated
% \ell_2 norm and the truncated \ell_1 norm.
% 
%   [y,final] = PROJRNORMBALL(z,r,p) produces the projection of z onto the 
%   unit-ball of the truncated \ell_p norm: ||.||_{\ell_p,r}, i.e. 
%   a) y is the solution to the optimization problem: 
%
%       minimize_X ||y-z||_{\ell_2} s.t. ||y||_{\ell_p,r} <= 1,
%
%   where p=1 or p=2.
%   b) final.t and final.s are the final values of parameters to the two 
%   nested search. For p = 1, final.k is the final value of the third inner 
%   search parameter.
%   If ||z||_{\ell_p,r} <= 1, then t = s = k = [].
%
%   [...] = PROJRNORMBALL(z,r,p,option) allows to specify furhter options:
%       1. [...] = PROJRNORMBALL(z,r,p,...,'tol',tol) sets the relative
%       tolerance of the deciding about zeros e.g.
%       For all i: y_i = 0 if |y_i| <= tol ||y||_{\ell_p,r). 
%       Default value: tol = 1e-12.
%       2. [...] = PROJRNORMBALL(z,r,p,...,'search',search) changes 
%       from default binary search to linear search over
%       a) t if search.t = 0.
%       b) s if search.s = 0.
%       c) k if search.k = 0 (only for p =1).
%       3. [...] = PROJRNORMBALL(z,r,p,...,'init',init) changes 
%       from default binary search start values t_0 = 1, k_0 = 1, s_0 = 0
%       to 
%       a) t_0 if init.t = t_0.
%       b) s_0 if init.s = s_0.
%       c) k_0 if init.k = k_0 (only for p =1).
 
function [y,final] = projrnormball(z,r,p,varargin)
n = length(z);
% Default values:
%
tol = 1e-12;  % relative tolerance
search.t = 1; % Flag binary search for t 
search.s = 1; % Flag binary search for s 
search.k = 1; % Flag binary search for k

% Set intial values for the searche parameters
init.t = 1; % Start value for t in the search
init.s = 0; % Start value for s in the search
init.k = 1; % Start value for k in the search
    

%% Check options
 
for i = 2:2:length(varargin)
     
    if strcmp(varargin{i-1},'tol')
        tol = varargin{i}; % Overwrite default tolerance
     
    elseif strcmp(varargin{i-1},'search') 
         
        % Flag linear search
        if isfield(search,'t') % Flag linear search for t
            search.t = varargin{i}.t;
        end
        if isfield(search,'s') % Flag linear search for s
            search.s = varargin{i}.s;
        end
        if p == 1 && isfield(search,'k') % Flag linear search for k if p = 1
            search.k = varargin{i}.k;
        elseif p ~= 1
            disp('k is only searched for if p = 1')
        end
    elseif strcmp(varargin{i-1},'init')
        % Set 
        if isfield(init,'t') % Define t_0
            if varargin{i}.t >= 1 && varargin{i}.t <= r 
                init.t = varargin{i}.t;
            else
                error('Choose t_0 such that 1<=t_0<=r!!!')
            end
        end  
        if isfield(init,'s') % Define s_0
            if varargin{i}.s >= 0 && varargin{i}.s <= n-r
                init.s = varargin{i}.s;
            else
                 error('Choose s_0 such that 0<=s_0<=length(z)-r!!!')
            end
         end
        if p == 1 && isfield(init,'k') % Define k_0 if p = 1
            if varargin{i}.k >= 1 && varargin{i}.k <= r
                init.k = varargin{i}.k;
            else
                error('Choose k_0 such that 1<=k_0<=r-t_0+1!!!')
            end
         elseif p ~= 1
            disp('k_0 can only be defined if p = 1')
        end
    end 
end  
  
% Intialize y
y = zeros(n,1);

% Store original sign of z.
I_sign = (z < 0);

% Replace z by sorted abs(z) and store orgininal sorting in Io.
[z,Io] = sort(abs(z),'descend');

% Pre-compute necessary sums dependend on t
if p == 2
    % Compute sum(z(r-t+1:r).^2) for all t
    sum_z_square = zeros(r,1);
    sum_z_square(1) = z(1)^2;
    for i=2:r
        sum_z_square(i) = sum_z_square(i-1)+z(i)^2;
    end
    tol_rel = sum_z_square(r)*tol; % Relative toleranz for the value 0
    if sum_z_square(r) - 1<= tol_rel 
        y = z;
        final.t = [];
        final.s = [];
        run_search = 0; %Flag: no search needed
    else
        run_search = 1; %Flag to search
    end
else
    % Compute sum(z(r-t+1:r)) for all t
    sum_z_t = zeros(r,1);
    sum_z_t(1) = z(r);
    for i=2:r
        sum_z_t(i) = sum_z_t(i-1)+z(r-i+1);
    end
    tol_rel = sum(z(1:r))*tol; % Relative toleranz for the value 0
    if sum_z_t(r)-1 <= tol_rel
        y = z;
        final.t = []; 
        final.s = [];
        final.k = [];
        run_search = 0; % Flag: no search needed
    else
        run_search = 1; % Flag to search
    end
end

if run_search == 1 
    
% Initialize binary search summation for computing z_{r-t+1} based on
% changes in s (and t for p = 2)
if search.t == 1 && p == 2
    sum_z_t = zeros(r,1);
    t_old = 1; % largest of all previous iterations
    sum_z_t(1) = z(r); 
end

if search.s == 1
    s_old = 1; % largest of all previous iterations
    sum_z_s = zeros(n-r+2,1); % sums for i=r+1,...,s
    if r+1<=n
        sum_z_s(2) = z(r+1); %Index shift because of s=0
    end
end

%% Search for t
t = init.t;
indt = 0;
 
t_min = 1; % Lower limit on t
t_max = r; % Upper limit on t

if search.t == 0
    init.t = 1;
    sum_zt = sum(z(r-init.t+1:r)); %sum_z_t for linear search
end
 
while indt == 0

    if t_min == t_max
        indt = 1;
        t = t_min;
    end
    
    if p == 2
        if search.t == 0
            if t_min ~= t_max && t > init.t%sum_zt was already update in the iteration before
                sum_zt = sum_zt + z(r-t+1);
            end
        else
            % Binary search summation for sum_z_t
            if t_old < t
                for i = t_old+1:t
                    sum_z_t(i) = sum_z_t(i-1) + z(r-i+1); 
                end
            t_old = t;
            end
            sum_zt = sum_z_t(t);
        end
    else
        sum_zt = sum_z_t(t);
    end
    
    %% Binary search for s
    s_min = 0;
    s_max = n-r;
    s = init.s;
    inds = 0;
    if search.s == 0
        init.s = 0;
        sum_zs = sum(z(r+1:r+init.s));
    end
    
        
  
    
    while inds == 0
        
        if s_min == s_max
            inds = 1;
            s = s_max;
        end
        
        %% Binary search summation for sum_z_s
        if search.s == 0
            if s_max ~= s_min && s>init.s %Only sums from r+1,...,s and sum_zs already update in the previous iteration
                sum_zs = sum_zs + z(r+s);
            end
        else
            if s_old < s
                for i = s_old+1:s
                    sum_z_s(i+1) = sum_z_s(i) + z(r+i); %Index shift by +1, because of s=0
                end
                s_old = s;
            end
            sum_zs = sum_z_s(s+1);
        end
        
        % Determine c_2 for p = 2 and z_tilde for both cases
        c_2 = (sum_zs + sum_zt);
        z_tilde = c_2/sqrt(t+s); %z_tilde_{r-t+1}
        
        % Determine mu    
        if p == 2
            
            c_1 = 0;
            if r-t > 0
                c_1 = sum_z_square(r-t);
            end
            % Use mu_shift = 1+mu. Then the polynomial in mu_shift is given by
            poly = [t^2, 2*s*t, -c_2^2*t+s^2-c_1*t^2, -2*c_1*s*t, -c_1*s^2];
            R = roots(poly);
            I = (abs(imag(R)) <= tol_rel);
            mu_shift = max(R(I)); % Choose largest real root
            mu = mu_shift-1; 
            
        else % Determine mu for p == 1
            
            % Sort the break points
            z_hat = zeros(r,1);
            z_hat(1:r-t) = z(1:r-t); %z_tilde_{1:r-t} = z_{1:r-t}
            z_hat(r-t+1) = t/sqrt(t+s)*z_tilde;
            [z_hat,I_alpha] = sort(z_hat,'descend');
        
            alpha = ones(r,1);
            alpha(r-t+1) = t^2/(t+s);
            alpha = alpha(I_alpha);
            
            % Binary search for k
            indk = 0;
        
            % Initialize binary search variables
            k_min = 1;
            k_max = r-t+1;
            k = init.k;
            k_old = 1;
             if search.k == 0
                init.k = 1;
                sum_z_hat = sum(z_hat(1:init.k));
                sum_alpha = sum(alpha(1:init.k));
             else
                % Intialize binary search summation of z_hat
                sum_z_hat_k = zeros(r,1);
                sum_z_hat_k(1) = z_hat(1);
        
                % Intialize binary search summation for alpha
                sum_alpha_k = zeros(r,1);
                sum_alpha_k(1) = alpha(1);
             end
        
            while indk == 0
                         
                if k_min == k_max
                    indk = 1;
                    k = k_max;
                end
            
                % Binary search summation for sum_z_hat_k and sum_alpha
                if search.k == 0
                   if k_max ~= k_min && k>init.k
                   sum_z_hat = sum_z_hat + z_hat(k);
                   sum_alpha = sum_alpha + alpha(k);
                   end
                else
                    if k_old < k   
                        for i=k_old+1:k
                            sum_z_hat_k(i) = sum_z_hat_k(i-1)+z_hat(i);
                            sum_alpha_k(i) = sum_alpha_k(i-1)+alpha(i);
                        end
                        k_old = k;
                    end
                    sum_alpha = sum_alpha_k(k);
                    sum_z_hat = sum_z_hat_k(k);
                end
                mu = (sum_z_hat-1)/sum_alpha;
            
                % Search rules for k
                if (z_hat(k) - alpha(k)*mu) >= -tol_rel
                    if mu >= -tol_rel
                        if search.k == 0
                            k_min = k;
                            k = k_min+1;
                        else
                            k_min = k;
                            k = k_min+ceil((k_max-k_min)/2); %ceil makes sure to not stay
                        end
                    else
                        k_min = k+1;
                        if search.k == 0
                            k = k_min;
                        else
                            k = k_min+floor((k_max-k_min)/2); %floor makes sure not to jump too far
                    
                        end
                    end
                else
                    if search.k == 0
                        k_min = k_max;
                    else
                        k_max = k-1; 
                        k = k_max-floor((k_max-k_min)/2); %floor makes sure to not jump too far, because k_max = k-1
                    end
                end
                      
            end
        end
        % Determine [y_{r+s}^(s,t),y_{r+s+1}^(s,t)] = [y_tilde_{r-t+1}/sqrt{s+t},z_{r+s+1}].
        if p == 2
            y(r+s) = (z_tilde/(1+mu*t/(s+t)))/sqrt(s+t);
        else
            y(r+s) = max(z_tilde-t*mu/sqrt(t+s),0)/sqrt(s+t);
        end
        if r+s < n
        y(r+s+1) = z(r+s+1);
        end
        
        % Search rules for s
        if r+s < n
            if y(r+s) - y(r+s+1) >= -tol_rel
                if search.s == 0
                    s_max = s_min;
                else
                    s_max = s;
                    s = s_max - ceil((s_max-s_min)/2); 
                end
            else
                s_min = s+1;
                if search.s == 0
                    s = s_min; 
                else
                    s = s_min + floor((s_max-s_min)/2);
                end
            end
        else
           if search.s == 0
               s_min = s_max;
           else
               s_max = s;
               s = s_max - ceil((s_max-s_min)/2);
           end
        end
    end
        % Determine [y_{r-t}^(s,t), y_{r-t+1}^(s,t)] = [y_tilde_{r-t},y_tilde_{r-t+1}/sqrt{s+t}]
      
        if r-t>0 % In case r = t there is nothing to do
            if p == 2
                y(r-t) = z(r-t)/(1+mu);
            else
                y(r-t) = max(z(r-t)-mu,0);
            end
        end
        y(r-t+1) = y(r+s);
        
        % Search rules for t
        if r-t > 0 
            
            if y(r-t) - y(r-t+1) >= -tol_rel
                if search.t ==0
                    t_max = t_min; %Tried all smaller t >= t is optimal
                else
                    t_max = t;
                    t = t_max - ceil((t_max-t_min)/2); %t may be too large
                end
            else 
                t_min = t+1;
                
                if search.t == 0
                    t = t_min; % Tried all previous ones => Only solution reached
                else
                    t = t_min + floor((t_max-t_min)/2); % t may be too large
                end
            end
        else
            if search.t == 0
                t_min = t_max;
            else
                t_max = t;
                t = t_max - ceil((t_max-t_min)/2); 
            end
        end
end
%% Get full solution
if p == 2
    y(1:r-t) = z(1:r-t)/(1+mu);
    y(r-t+1:r+s) = (z_tilde/(1+mu*t/(s+t)))/sqrt(s+t);
    y(r+s+1:end) = z(r+s+1:end);
else
    y(1:r-t) = max(z(1:r-t)-mu,0);
    y(r-t+1:r+s) = max(z_tilde-t*mu/sqrt(t+s),0)/sqrt(s+t);
    y(r+s+1:end) = z(r+s+1:end);
end

    final.t = t;
    final.s = s;
    if p == 1
        final.k = k;
    end
end
y(y<tol_rel) = 0;


[~,Iback] = sort(Io);
y = y(Iback);

y(I_sign) = -y(I_sign);
end