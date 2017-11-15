function reprojErr = computeReprojectionError(Ps, X, images)
% COMPUTEREPROJECTIONERROR  Function that takes in a set of projection
% matrices, 3D, and 2D points, and computes the reprojection error. Assumes 
% a general 3 x 4 projection matrix and computes reprojection error as the
% expression (P*X - x), where P is the projection matrix (3 x 4), X is a
% homogeneous 3D point (4 x 1), and x is a homogeneous image point (3 x 1).
% The function returns a reprojection error vector of an appropriate size,
% which contains reprojection errors of the x and y image coordinates for
% each 3D point in each image stacked as follows. If each image contains N
% points, and there are M images, then the reprojection error vector is of
% size 2*M*N x 1, and it first contains all 2*N errors (x and y image
% coordinates) for image 1, all 2*N errors for image 2, and so on up to
% image M.


% Number of views
numViews = size(Ps,1);
% Number of (3D) points
numPoints = size(X,2);

% Initialize the reprojection error vector
reprojErr = zeros(2 * numViews * numPoints, 1);

% Compute reprojection errors
for i = 1:numViews
    % Project 3D points to 2D
    x_hat = squeeze(Ps(i,:,:)) * X;
    x_hat = x_hat ./ repmat(x_hat(3,:), 3, 1);
    % Compute repojection error with respect to the observed 2D points
    errTemp = x_hat(1:2,:) - images{i}(1:2,:);
    errTemp = errTemp(:);
    % Set small errors to zeros
    errTemp(find(abs(errTemp) <= 1e-5)) = 0.0;
    reprojErr((i-1)*2*numPoints+1 : 2*i*numPoints) = errTemp;
end

end
