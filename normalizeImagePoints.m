function [normalizedPoints, normalizingTransform] = normalizeImagePoints(unnormalizedPoints)
% NORMALIZEIMAGEPOINTS  Function to normalize image points such that the
% mean of the points is the origin and that the mean of the distance of the
% normalized points from the origin is sqrt(2). This kind of normalization
% is widely used to avoid numerical instability of linear systems.
% The function also returns the normalizing transform, so that it can be
% used for 'unnormalization'.


% Initialize the variable that stores the result
normalizedPoints = unnormalizedPoints;

% Center data
centroid = mean(unnormalizedPoints, 2);
normalizedPoints(1:2,:) = unnormalizedPoints(1:2,:) - repmat([centroid(1); centroid(2)], 1, size(unnormalizedPoints,2));

% Compute the mean distance of points from the origin
d = sum(sqrt(sum(normalizedPoints(1:2,:).^2))) / size(normalizedPoints,2);

% If the average distance is zero, return the current set of points
if d == 0
    fprintf('Warning: The supplied points are zeros.\n');
    return;
end

% Compute the scaling factor
s = sqrt(2) / d;

% Use the scaling factor to scale the points
scalingTransform = [s, 0, 0; 0, s, 0; 0, 0, 1];
normalizedPoints = scalingTransform * normalizedPoints;

% Also compute the normalizing transform
normalizingTransform = [s, 0, -s*centroid(1); 0, s, -s*centroid(2); 0, 0, 1];

end
