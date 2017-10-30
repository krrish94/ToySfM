%% Script to test the ToySfM example


%% Synthesize the scene
% The scene comprises a cube (whose center also happens to be the world
% origin). A camera travels around the cube, capturing (usually 8) images
% of the cube from diverse viewpoints. We assume some (perfectly) known
% camera intrinsics and use that to synthesize images captured by the
% camera.

% Number of views (= num of cameras = num of images)
numViews = 8;

% Camera intrinsics (focal length, retina-to-image coordinate shifts, focal
% length scaling). Assume fx = fy.
f = 1;
alpha_ccd = 100;
cx = 50;
cy = 100;
% Intrinsics matrix
K = [alpha_ccd*f, 0, cx; 0, alpha_ccd*f, cy; 0, 0, 1];

% Generate the cube (centered at the origin; length of each side is 4
% units; each edge contains 4 points)
% Ground-Truth world points
worldPoints_gt = generateCube(4,4);

% Simulate the camera trajectory (circle) around the cube
[Rs, ts] = generateCameraTrajectory(numViews, 10);

% Camera locations
camLocs = ts';

% % Visualize the scene and the camera trajectories
% scatter3(worldPoints_gt(1,:), worldPoints_gt(2,:), worldPoints_gt(3,:), 'filled');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis('equal');
% % axis([-5,5,-5,5,-5,5]);
% hold on;
% scatter3(camLocs(1,:), camLocs(2,:), camLocs(3,:), 'filled');

% Synthesize images of the cube
images = synthesizeImages(worldPoints_gt, K, Rs, ts);

% % Visualize images and camera coordinates
% for i = 1:numViews
%     scatter(images{i}(1,:), images{i}(2,:), 'filled');
%     xlim([0 100]);
%     ylim([50 150]);
%     xlabel('y');
%     ylabel('x');
%     pause(0.5);
% end


%% Construct the initial guess

% Use the first two views to compute the Fundamental matrix

