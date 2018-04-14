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
function X = projindex(Z,N,Index,H,gamma)
X = Z - gamma*H;
X(Index) = N(Index);    
end