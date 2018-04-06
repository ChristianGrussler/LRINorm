%%%% Example: Model Reduction of a Linear Discrete Time SISO System
% Requirements to run this file: CVX (and SDP solver)
% 
%   Steps: 
%   1. Approximate finite dimensional Hankel matrix of the system with 
%      low-rank inducing Frobenius norm (r* norm) in CVX
%   2. Obtain realization that matches the impulse response of the
%      approximation with Kung's realization algorithm
%
%%%%

% Choose desired order r
r = 1;

% Define SISO system matrices
A = diag(0:.1:.9);
n = length(A);
B = ones(n,1);
C = ones(n,1)';
D = 0;

% Define state-space system
G = ss(A,B,C,D,1);

%% r* approach
% Determine impulse response 
h = impulse(G,14*n+1);
n = length(h);

% Transform into Hankel matrix (without t = 0)
H = hankel(h(2:n/2+1),h(n/2+1:n));

% Define SDP program in CVX for r*norm approximation
[n,m] = size(H);
cvx_begin SDP 
cvx_precision best
            variable W(n,n) semidefinite 
            variable P(m,m) semidefinite
            variable M(n,m) hankel
                minimize 1/2*trace(P)-trace(H'*M)
                subject to
                [eye(n) - W M; transpose(M) P] >= 0;
                W >= 0;
                trace(W) == n - r;
cvx_end

% Extract impulse response of reduced system from solution M
hr = [G.D; M(:,1); M(end,2:end)'];

% Apply Kung's realization algorithm to hr
n = length(hr);
l = floor((n-1)/2);
Ho = hankel(hr(2:l+2),hr(l+2:2*l+1));
Hb = hankel(hr(3:l+3),hr(l+3:2*l+2));

[U,S,V] = svd(Ho);
S = diag(S);
Cr = sqrt(diag(S(1:r)))*V(:,1:r)';
Or = U(:,1:r)*sqrt(diag(S(1:r)));
Ghr = ss(Or\Hb/Cr,Cr(:,1),Or(1,:),hr(1),G.ts);

%% Determine balanced truncated model Gr of order r
Gb = balreal(G);
A = Gb.a;
B = Gb.b;
C = Gb.c;
Gr = ss(A(1:r,1:r),B(1:r,:),C(:,1:r),G.D,G.ts); % Note that balred gives sometimes worse approximations!

%% Relative Error comparison of the H_infty norm
e_rstar_vs_bal = norm(G-Ghr,'inf')/norm(G-Gr,'inf');