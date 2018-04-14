% Example: Matrix completion by low-rank Frobenius-norm mimization
% N = svd_5(H); 
% Known entries: positive elements
% Unkown entries: non-positive elements
% Compare results of DR, NDR and CVX to illustrate their handlings
%
%%%%%%%%%%%%%
% References: 
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
clc
H = hankel(ones(10,1));
r = 5;
[U,S,V] = svd(H);
N = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
[n,m] = size(N);
% Define known entries
Index = H > 0;

% Intialize Douglas-Rachford iterations in the same point: Z0 = 0
Z0 = zeros(n,m);
% Set iteration tolerance smaller than default
tol = 1e-9;

%% Compute the different solutions:

% Douglas-Rachford solving the convexified problem
tic;
[M_dr,rankM_dr,err_dr,D_dr,Z_fix_dr,iter_dr] = drcomplete(N,Index,r,2,'Z0',Z0,'tol',tol);
toc; t_dr = toc;

% Non-convex Douglas-Rachford
tic;
[M_ndr,rankM_ndr,err_ndr,D_ndr,Z_fix_ndr,iter_ndr] = drcomplete(N,Index,r,2,'solver','NDR','Z0',Z0,'tol',tol);
toc; t_ndr = toc;

% CVX
tic;
[M_cvx,rankM_cvx,err_cvx,D_cvx] = cvxcomplete(N,Index,r,2);
toc; t_cvx = toc;

%% Display summary
disp('Rank of the different solvers:');
disp('CVX');
disp(rankM_cvx);
disp('Dougals-Rachford');
disp(rankM_dr);
disp('Non-convex Dougals-Rachford');
disp(rankM_ndr);

pause

disp('Relative Erros of the different solvers:');
disp('CVX');
disp(err_cvx/norm(N,'fro'));
disp('Dougals-Rachford');
disp(err_dr/norm(N,'fro'));
disp('Non-convex Dougals-Rachford');
disp(err_ndr/norm(N,'fro'));

pause 

disp('Elapse time of different solvers:');
disp('Dougals-Rachford');
disp(t_dr);
disp('Non-convex Dougals-Rachford');
disp(t_ndr);
disp('CVX');
disp(t_cvx);
