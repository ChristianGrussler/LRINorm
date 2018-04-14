%PROJINDEX Orthogonal projection onto the subspace of known entries
%
%   X = PROJHANKEL(Z,N,Index,H,gamma) determines the orthogonal projection of 
%   Z-gamma*H onto the subspaces of matrices with entries N(Index), i.e., 
%   X is the proximal mapping of the function
% 
%   gamma*(i_{X(Index) = N(Index)}(X)+trace(X'H)),
%
%   where i_{X(Index) = N(Index)} is the indicator function for the linear
%   constraint X(Index) = N(Index).   
%
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
function X = projindex(Z,N,Index,H,gamma)
X = Z - gamma*H;
X(Index) = N(Index);    
end