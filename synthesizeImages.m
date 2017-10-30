function images = synthesizeImages(X, K, Rs, ts)
% SYNTHESIZEIMAGES  Takes in a 3 x N vector of N 3D points (X), a 3 x 3
% camera intrinsics matrix (K), and a set of M rotation matrices (M x 3 x
% 3) (Rs), and corresponding translation vectors (M x 3) (ts). Synthesizes
% M images, such that each image contains projections of all the N 3D
% points. The synthesized images are stored in a cell array of length M,
% and each entry in the cell arry is a 2 x N matrix containing the x and y
% image coordinates of the projections of each 3D point. While converting
% from world to camera coordinates, we assume that the convention is
% X_cam = R^T * (X_world - t)


% Number of views
numViews = size(ts,1);

% Number of points
numPoints = size(X,2);

% Cell array to store images
images = cell(numViews);

% Synthesize the images
for i = 1:numViews
    images_homogeneous = K * (squeeze(Rs(i,:,:))' * (X - repmat(ts(i,:)',1,numPoints)));
    images_homogeneous = images_homogeneous ./ repmat(images_homogeneous(3,:), 3, 1);
    images{i} = images_homogeneous(1:2,:);
end


end