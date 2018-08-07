% Example: Low-rank approximation with Hankel constraint
% Compare results of DR, NDR and CVX to illustrate their handlings
%
%%%%%%%%%%%%%
% References: 
%   - C. Grussler and P. Giselsson (2017):
%   "Local convergence of proximal splitting methods for rank constrained
%   problems", pp. 702-708, IEEE 56th Annual Conference on Decision and Control
%   (CDC), DOI: 10.1109/CDC.2017.8263743.
%
%   - C. Grussler and P. Giselsson (2016):
%   "Low-Rank Inducing Norms With Optimality Interpreations", 
%   arXiv:1612.03186v1.
%
%   - C. Grussler (2017):
%   "Rank reduction with convex constraints", PhD Thesis, 
%   Department of Automatic Control, Lund Institute of Technology, 
%   Lund University, ISBN 978-91-7753-081-7.  
%%%%%%%%%%%%%

clc

H = hankel(1:10,10:-1:1); % Hankel matrix
r = 5; % Desired rank of the approximation
[n,m] = size(H);
% Compute the different solutions:
% Set tolerance for deciding about multiple singular values
tol = 1e-10;

% Low-rank inducing Frobenius norm with CVX
tic;
[M_cvx,rankM_cvx,err_cvx,D_cvx] = cvxhankelapprox(H,r,'tol',tol);
toc; t_cvx = toc;

% Low-rank inducing Frobenius norm with Douglas-Rachford
tic;
[M_dr,rankM_dr,err_dr,D_dr,Z_fix_dr,iter_dr] = drhankelapprox(H,r,'tol',tol);
toc; t_dr = toc;

% Non-convex Douglas-Rachford
tic;
[M_ndr,rankM_ndr,err_ndr,D_ndr,Z_fix_ndr,iter_ndr] = drhankelapprox(H,r,'solver','NDR','tol',tol);
toc; t_ndr = toc;

% Display summary
disp('Rank of the solutions:');
disp('CVX');
disp(rankM_cvx);
disp('Dougals-Rachford');
disp(rankM_dr);
disp('Non-convex Dougals-Rachford');
disp(rankM_ndr);

pause

disp('Relative Errors of the solutions:');
disp('CVX');
disp(err_cvx/norm(H,'fro'));
disp('Dougals-Rachford');
disp(err_dr/norm(H,'fro'));
disp('Non-convex Dougals-Rachford');
disp(err_ndr/norm(H,'fro'));

pause

disp('Elapse time of the solvers:');
disp('Dougals-Rachford');
disp(t_dr);
disp('Non-convex Dougals-Rachford');
disp(t_ndr);
disp('CVX');
disp(t_cvx);
