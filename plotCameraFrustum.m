function [] = plotCameraFrustum(R, camLoc)
% PLOTCAMERAFRUSTUM  Plots the frustum of a camera that is rotated by R and
% translated by t with respect to the world coordinate frame.
% R describes the rotation matrix that takes the world coordinate frame to
% the camera coordinate frame. camLoc denotes the location of the camera
% centre in the world frame.

% Size unit of the frustum
frustumSize = 1;

% Vertices of the frustum. The first vertex is the camera center. Vertices
% 2, 3, 4, 5 are the top-left, top-right, bottom-right, and bottom-left
% corners of the image respectively.
frustumVertices = [0, 0, 0;   -frustumSize, -0.5*frustumSize, frustumSize; ...
    frustumSize, -0.5*frustumSize, frustumSize;   ...
    frustumSize, 0.5*frustumSize, frustumSize;   ...
    -1*frustumSize, 0.5*frustumSize, frustumSize]';

% Convert the frustum to world coordinates
frustumVertices = R*frustumVertices + repmat(camLoc, 1, 5);

% Plot the frustum
frustumEdges = [1,2; 1,3; 1,4; 1,5; 2,3; 3,4; 4,5; 5,2];
for i = 1:size(frustumEdges,1)
    plot3(frustumVertices(1,[frustumEdges(i,1), frustumEdges(i,2)]), ... 
        frustumVertices(2,[frustumEdges(i,1), frustumEdges(i,2)]), ...
        frustumVertices(3,([frustumEdges(i,1), frustumEdges(i,2)])), 'b');
end

end
