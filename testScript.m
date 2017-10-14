%% Script to test the ToySfM example


%% Synthesize the scene
% The scene comprises a cube (whose center also happens to be the world
% origin). A camera travels around the cube, capturing (usually 8) images
% of the cube from diverse viewpoints. We assume some (perfectly) known
% camera intrinsics and use that to synthesize images captured by the
% camera.

% Generate the cube (centered at the origin; length of each side is 4
% units; each edge contains 4 points)
% Ground-Truth world points
worldPoints_gt = generateCube(4,4);

% Simulate the camera trajectory (circle) around the cube
[Rs, ts] = generateCameraTrajectory(8, 10);

% Visualize the scene and the camera trajectories
scatter3(worldPoints_gt(1,:), worldPoints_gt(2,:), worldPoints_gt(3,:), 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis([-5,5,-5,5,-5,5]);
