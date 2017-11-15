function [F,A] = computeFundamentalMatrix(matches1, matches2)
% COMPUTEFUNDAMENTALMATRIX  Takes as input two 3-by-N matrices containing
% point correspondences between two images (in homogeneous coordinates) and
% computes the fundamental matrix F.


% Normalize image points
% [img1, normTf1] = normalizeImagePoints(matches1);
% [img2, normTf2] = normalizeImagePoints(matches2);
img1 = matches1;
img2 = matches2;

% Form the linear system Af = 0
% A = [img2(1,:)'.*img1(1,:)', img2(1,:)'.*img1(2,:)', img2(1,:)'.*img1(3,:)', ...
%     img2(2,:)'.*img1(1,:)', img2(2,:)'.*img1(2,:)', img2(2,:)'.*img1(3,:)', ...
%     img1(1,:)', img1(2,:)', img1(3,:)'];
A = [img2(1,:)'.*img1(1,:)', img2(1,:)'.*img1(2,:)', img2(1,:)', ...
    img2(2,:)'.*img1(1,:)', img2(2,:)'.*img1(2,:)', img2(2,:)', ...
    img1(1,:)', img1(2,:)', ones(size(img1,2),1)];


% Solve the system using SVD and take f to be the eigenvector corresponding
% to the smallest eigenvalue, i.e., take f to be the last column of V.
[~, ~, V] = svd(A, 0);
f = V(:,9);
F_ = reshape(f, 3, 3)';
clear f;

% Enforce the rank-2 constraint explicitly
[U, D, V] = svd(F_, 0);
D
F = F_;
% F = U * diag([D(1,1), D(2,2), 0]) * V';
clear F_ U D V;

% % Unnormalize
% F = normTf2' * F * normTf1;

end
