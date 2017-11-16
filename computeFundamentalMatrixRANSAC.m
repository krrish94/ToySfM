function F = computeFundamentalMatrixRANSAC(matches1, matches2, maxIters)
% COMPUTEFUNDAMENTALMATRIXRANSAC  Takes as input two 3-by-N matrices 
% containing point correspondences between two images (in homogeneous 
% coordinates) and computes the fundamental matrix F. Uses a RANSAC routine
% to discard outliers


% If maxIters is not specified, set it to 50
if nargin < 3
    maxIters = 50;
end

% Largest Inlier Set
inliers = [];
numInliers = -Inf;

% Inlier threshold
inlierThresh = 0.01;

% Main RANSAC loop
for i = 1:maxIters
    
    % Randomly sample 8 points
    inds = randsample(size(matches1,2),8);
    
    % Use these samples to compute the fundamental matrix
    F_cur = computeFundamentalMatrix(matches1(:,inds), matches2(:,inds));
    
    % Inlier set (current)
    inliers_cur = [];
    numInliers_cur = 0;
    
    % Estimate the number of inliers
    for j = 1:size(matches1,2)
        if abs(matches2(:,j)' * F_cur * matches1(:,j)) <= inlierThresh
            inliers_cur = [inliers_cur, j];
            numInliers_cur = numInliers_cur + 1;
        end
    end
    
    % Best inilier set so far
    if numInliers_cur > numInliers
        numInliers = numInliers_cur;
        inliers = inliers_cur;
    end
    
    % If more than 50 percent of the data are inliers, we assume a
    % reasonable estimate has been found and terminate.
    if numInliers >= 0.5 * size(matches1,2)
        break
    end
    
end

% Refine the fundamental matrix estimate using the largest inlier set
F = computeFundamentalMatrix(matches1(:,inliers), matches2(:,inliers));

end
