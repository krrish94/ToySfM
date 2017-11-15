function [R, t] = computeRelativePose(image1, image2, K)
% COMPUTERELATIVEPOSE  Takes in two images (homogeneous image coordinates)
% and computes the relative pose between the two images. The computed
% relative pose R, t, when converted to a homogeneous 4 x 4 transform
% matrix as [R, t; zeros(1,3), 1], transforms 3D points from the camera
% coordinates of the camera capturing the first image to the camera
% coordinates of the second image. Assumes that the intrinsics (K) are 
% known.


% Use the first two views to compute the fundamental matrix. The function
% 'computeFundamentalMatrix' is used for this purpose. This function
% expects homogenized image coordinates, hence we append a 1 to the image
% coordinates before passing them to this function.
F_12 = computeFundamentalMatrix(image1, image2);
% F_12 = fundmatrix(image1, image2);
% F_12 = optimalRansacfitfundmatrix(image1, image2, 50);
% F_12 = computeFundamentalMatrixRANSAC(image1, image2);
% Get the essential matrix from the fundamental matrix
E_12 = K' * F_12 * K;
% Get the closest rank-2 approximation to the essential matrix
[U_12, D_12, V_12] = svd(E_12, 0);
E_12 = U_12 * diag([D_12(1,1), D_12(1,1), 0]) * V_12;
% Clear temporary variables
clear U_12 D_12 V_12

% Decompose the essential matrix to obtain relative pose (rotation,
% translation) between the two views. The rotation (R) and translation (t)
% thus obtained are in a way that the center of camera 1 (C1) can be
% expressed with respect to the center of camera 2 (C2) as per the equation
% C2 = R * C1 + t.
[R, t] = decomposeEssentialMatrix(E_12, image1, image2, K);

end
