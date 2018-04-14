%PROJHANKEL Orthogonal projection onto the subspace of Hankel matrices
%
%   X = PROJHANKEL(Z,H,gamma) determines the orthogonal projection of 
%   Z-gamma*H onto the subspaces of Hankel matrices, i.e., X is the 
%   proximal mapping of gamma*(i_{Hankel}(X)+trace(X'H)), where i_{Hankel} 
%   is the indicator function of the set of Hankel matrices.
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
function X = projhankel(Z,H,gamma)

Z = Z - gamma*H;
dim = size(Z);
if dim(1) < dim(2)
    Z = Z';
    transp_Z = 1; % Flag transpositon of Z
    dim = sort(dim,'descend');
else
    transp_Z = 0; % Flag Z is not transposed
end
n = dim(1);
m = dim(2);

% Determine the means k(1),...,k(n+m-1) along the anti-diagonals
B = [Z;zeros(m)];
k = sum(reshape(B(1:(n+m-1)*m),n+m-1,m),2)./[1:m-1 m*ones(1,n-m+1) m-1:-1:1]'; 
% Arrange means as Hankel matrix
X = hankel(k(1:n),k(n:end));

if transp_Z == 1
    X = X';
end
    
    
end