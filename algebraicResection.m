function [P] = algebraicResection(X, x)
% ALGEBRAICRESECTION  Takes in a set of 3D points (4 x N matrix) and a set
% of homogeneous image points (3 x N matrix), and computes the projection
% matrix P (3 x 4 matrix) using a Direct Linear Transformation (DLT)
% method. First, we set up a least squares system, and then solve it by
% computing the SVD of the coefficient matrix and picking the eigenvector
% corresponding to the smallest non-zero eigenvalue as our solution.


% Coefficient matrix
A = zeros(2*size(X,2), 12);

% Construct the linear system
for i = 1:size(X,2)
    A(2*i-1,:) = [X(:,i)'*x(3,i), zeros(1,4), -1*X(:,i)'*x(1,i)];
    A(2*i,:) = [zeros(1,4), X(:,i)'*x(3,i), -1*X(:,i)'*x(2,i)];
end

% Solve, using SVD
[~, D, V] = svd(A,0);
p = V(:,end-2);
p = p ./ p(end);
P = reshape(p,[4,3])';

end
