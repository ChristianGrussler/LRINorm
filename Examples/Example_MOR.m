%%%% Example: Model Reduction of a Linear Discrete Time SISO System
% 
%   Steps: 
%   1. Approximate finite dimensional Hankel matrix of the system with 
%      low-rank inducing Frobenius norm (r* norm) using
%      a) CVX (Requires CVX and SDP solver!).
%      b) Douglas-Rachford 
%      and with (non-convex) Douglas-Rachford.
%   2. Obtain realization that matches the impulse response of the
%      approximation with Kung's realization algorithm
%   
%   Comparsion: 
%   i.)  (H-infty error of each method)/(H-infty error of Balanced
%   Truncation)
%   ii.) Elapse time of each method
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

% Choose desired order r
r = 5;

% Define SISO system matrices
A = diag(0:.1:.9);
n = length(A);
B = ones(n,1);
C = ones(n,1)';
D = 0;

% Define state-space system
G = ss(A,B,C,D,1);


% Determine impulse response 
h = impulse(G,14*n+1);
n = length(h);

% Transform into Hankel matrix (without t = 0)
H = hankel(h(2:n/2+1),h(n/2+1:n));


%% Approximate H with low-rank Hankel matrix:

% Convexification with low-rank inducing Frobenius norm using
% Douglas-Rachford:
tic;
[M_dr,rankM_dr,err_dr,D_dr] = drhankelapprox(H,r,'tol',1e-12); 
toc; 
t_dr = toc; % Elapse time for DR

% Extract impulse response from approximations:
hr_dr = [G.D; M_dr(:,1); M_dr(end,2:end)'];

%% Apply Kung's realization algorithm to hr
G_dr = kung(hr_dr,r,G.ts);

%% Determine balanced truncated model Gr of order r
Gb = balreal(G);
A = Gb.a;
B = Gb.b;
C = Gb.c;
Gr = ss(A(1:r,1:r),B(1:r,:),C(:,1:r),G.D,G.ts); % Note that balred gives sometimes worse approximations!


%% Relative Error comparison of the H_infty norm
e_dr_vs_bal = norm(G-G_dr,'inf')/norm(G-Gr,'inf');
msg_dr = sprintf('(H-infty error of r* relaxation with DR)/(H-infty error BT): %s',e_dr_vs_bal);
disp(msg_dr);

pause
%% Approximate H with low-rank Hankel matrix:

% Use non-convex Douglas-Rachford:
tic;
[M_ndr,rankM_ndr,err_ndr,D_ndr] = drhankelapprox(H,r,'solver','NDR','tol',1e-12); 
toc; 
t_ndr = toc; % Elapse time for NDR

% Extract impulse response from approximations:
hr_ndr = [G.D; M_ndr(:,1); M_ndr(end,2:end)'];

%% Apply Kung's realization algorithm to hr
G_ndr = kung(hr_ndr,r,G.ts);


%% Relative Error comparison of the H_infty norm
e_ndr_vs_bal = norm(G-G_ndr,'inf')/norm(G-Gr,'inf');
msg_ndr = sprintf('(H-infty error of non-convex DR)/(H-infty error BT): %s',e_ndr_vs_bal);
disp(msg_ndr);

pause

%% Approximate H with low-rank Hankel matrix:
tic;
% Convexification with low-rank inducing Frobenius norm using CVX:
[M_cvx,rankM_cvx,err_cvx,D_cvx] = cvxhankelapprox(H,r); 
toc;
t_cvx = toc;

% Extract impulse response from approximations:
hr_cvx = [G.D; M_cvx(:,1); M_cvx(end,2:end)'];

%% Apply Kung's realization algorithm to hr
G_cvx = kung(hr_cvx,r,G.ts);

%% Relative Error comparison of the H_infty norm
e_cvx_vs_bal = norm(G-G_cvx,'inf')/norm(G-Gr,'inf');
msg_cvx = sprintf('(H-infty error of r* relaxation with CVX)/(H-infty error BT): %s',e_cvx_vs_bal);
disp(msg_cvx);

pause

%% Display Summary

disp('Elapse time of different solvers:');
disp('Dougals-Rachford');
disp(t_dr);
disp('Non-convex Dougals-Rachford');
disp(t_ndr);
disp('CVX');
disp(t_cvx);

pause 
disp('H-infty error/H-infty error of BT:');
disp('Dougals-Rachford');
disp(e_dr_vs_bal);
disp('Non-convex Dougals-Rachford');
disp(e_ndr_vs_bal);
disp('CVX');
disp(e_cvx_vs_bal);

