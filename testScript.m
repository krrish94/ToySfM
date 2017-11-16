%% Script to test the ToySfM example


%% Synthesize the scene
% The scene comprises a cube (whose center also happens to be the world
% origin). A camera travels around the cube, capturing (usually 8) images
% of the cube from diverse viewpoints. We assume some (perfectly) known
% camera intrinsics and use that to synthesize images captured by the
% camera.


% Seed the random number generator (for repeatability)
rng('default');

% Camera intrinsics (focal length, retina-to-image coordinate shifts, focal
% length scaling). Assume fx = fy (= f).
f = 1;
alpha_ccd = 10;
cx = 5;
cy = 10;
% Intrinsics matrix
% K = [alpha_ccd*f, 0, cx; 0, alpha_ccd*f, cy; 0, 0, 1];
K = eye(3);

% Generate the cube (centered at the origin; length of each side is 4
% units; each edge contains 4 points)
% Ground-Truth world points
worldPoints_gt = generateCube(4,16);
% worldPoints_gt = load('cubePts.mat', 'cubePts');
% worldPoints_gt = worldPoints_gt.cubePts;

% Number of views (= num of cameras = num of images)
numViews = 8;

% Number of points ( = number of points in the generated world)
numPoints = size(worldPoints_gt,2);

% Simulate the camera trajectory (circle) around the cube
[Rs, ts] = generateCameraTrajectory(numViews, 25);

% Camera locations
camLocs = ts';

% Initialize a variable to hold all projection matrices
Ps = zeros(numViews, 3, 4);

% Visualize the scene and the camera trajectories
scatter3(worldPoints_gt(1,:), worldPoints_gt(2,:), worldPoints_gt(3,:), 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('The synthetic scene and cameras');
axis('equal');
hold on;
for i = 1:numViews
    R_cam_W = squeeze(Rs(i,:,:));
    % Plot the camera frustum
    plotCameraFrustum(squeeze(Rs(i,:,:)), camLocs(:,i));
    % Plot the camera coordinate frame
    scatter3(camLocs(1,i), camLocs(2,i), camLocs(3,i), 25, [0, 0, 1], 'filled');
    text(camLocs(1,i), camLocs(2,i), camLocs(3,i)-2, sprintf('C%d',i), 'Color', 'red');
    quiver3(camLocs(1,i), camLocs(2,i), camLocs(3,i), R_cam_W(1,1), R_cam_W(2,1), R_cam_W(3,1), 'r', 'LineWidth', 2);
    text(camLocs(1,i) + R_cam_W(1,1), camLocs(2,i) + R_cam_W(2,1), camLocs(3,i) + R_cam_W(3,1), 'X', 'Color', 'red');
    quiver3(camLocs(1,i), camLocs(2,i), camLocs(3,i), R_cam_W(1,2), R_cam_W(2,2), R_cam_W(3,2), 'g', 'LineWidth', 2);
    text(camLocs(1,i) + R_cam_W(1,2), camLocs(2,i) + R_cam_W(2,2), camLocs(3,i) + R_cam_W(3,2), 'Y', 'Color', 'green');
    quiver3(camLocs(1,i), camLocs(2,i), camLocs(3,i), R_cam_W(1,3), R_cam_W(2,3), R_cam_W(3,3), 'b', 'LineWidth', 2);
    text(camLocs(1,i) + R_cam_W(1,3), camLocs(2,i) + R_cam_W(2,3), camLocs(3,i) + R_cam_W(3,3), 'Z', 'Color', 'blue');
    clear R_cam_W;
end



% Synthesize images of the cube
images = synthesizeImages(worldPoints_gt, K, Rs, ts);

% % Visualize the cube in camera coordinates
% figure;
% for i = 1:numViews
%     R = squeeze(Rs(i,:,:));
%     t = ts(i,:)';
%     X_cam = R' * (worldPoints_gt - t);
%     scatter3(X_cam(1,:), X_cam(2,:), X_cam(3,:), 'filled');
%     view(0,-90);
%     axis('equal');
%     pause(0.5);
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     title('Camera coordinates');
% end

% Visualize the formed images
% figure;
% for i = 1:numViews
%     scatter(images{i}(1,:), images{i}(2,:), 'filled');
%     % xlim([0 100]);
%     % ylim([50 150]);
%     xlabel('x');
%     ylabel('y');
%     title('Images formed');
%     pause(0.5);
% end


%% Construct the initial guess


% For each view (other than the first view, compute the relative pose of
% the camera with respect to the first view
R_cur = eye(3);
t_cur = zeros(3,1);
for i = 2:numViews
    
    % Compute the relatvie pose of the current view with respect to the
    % first view
    [R_21, t_21] = computeRelativePose(images{1}, images{i}, K);
    % [R_21, t_21, ~, ~, ~] = compute_relative_transformation(images{1}(1:2,:), images{i}(1:2,:), K);
    
    % If this is the second view, use this view to initialize the initial
    % structure estimate as well, by triangulation.
    % Triangulate using the first two views and their relative pose to obtain a
    % 3D structure initialization
    if i == 2
        
        Ps(1,:,:) = K * [eye(3), [0; 0; 0]];
        Ps(2,:,:) = K * [R_21, t_21];
        X_init = algebraicTriangulation(images{1}, images{2}, squeeze(Ps(1,:,:)), squeeze(Ps(2,:,:)));
        X_init = X_init ./ repmat(X_init(4,:), 4, 1);
        
        % % Visualize initial guess (triangulated from two views)
        % figure;
        % scatter3(X_init(1,:), X_init(2,:), X_init(3,:), 'filled');
        % xlabel('X');
        % ylabel('Y');
        % zlabel('Z');
        % title('Triangulated 3D structure using the first two views');
       
    % For all other views
    else
        
        % Construct the rotation matrix and translation vectors that
        % transform a point from the first frame (cam 1) to cam i frame
        % R_cur = R_cur * R_21;
        % t_cur = R_cur * t_21 + t_cur;
        R_cur = R_21;
        t_cur = t_21;
        
        % This translation has to be scaled by a scale factor. That scale
        % factor is computed as a ratio of reconstruction obtained by
        % triangulation of the first view and the current view and the
        % triangulation of the first view and the second view. For
        % instance, if a distance is measured as d1 in the triangulation of
        % the first two views, and as d2 in the triangulation of the first
        % view and the current view, the translation t_cur (i.e., t_21
        % computed above) should be rescaled as (d1/di)*t_21.
        P_temp = K * [R_21, t_21];
        X_temp = algebraicTriangulation(images{1}, images{i}, squeeze(Ps(1,:,:)), P_temp);
        X_temp = X_temp ./ repmat(X_temp(4,:), 4, 1);
        
        numSamplesForScale = 25;
        scaleFactors = zeros(numSamplesForScale,1);
        for j = 1:numSamplesForScale
            % Compute distance between two (arbitrary) points in the initial
            % triangulation (view 1 and view 2)
            sampledPoints = randsample(numPoints,2);
            d1 = sqrt(sum((X_init(1:3,sampledPoints(1)) - X_init(1:3,sampledPoints(2))).^2));
            % Compute distance between the same pair of points in the current
            % triangulation (view 1 and view i)
            di = sqrt(sum((X_temp(1:3,sampledPoints(1)) - X_temp(1:3,sampledPoints(2))).^2));
            % Compute the scaling factor
            scaleFactor = d1 / di;
            % Store the current scaling factor
            scaleFactors(j) = scaleFactor;
        end
        
        % Average of all scaling factors
        avgScaleFactor = mean(scaleFactors);
        
        % Scale the translation vector
        t_cur = avgScaleFactor * t_cur;
        
        % Construct the projection matrix using the estimated relative pose
        Ps(i,:,:) = K * [R_cur, t_cur];
        
    end
    
    % Clear temporary variables
    % clear R_21 t_21
    
end


% % % Perform resection to recover the projection matrices for the remaining
% % % views
% % for i = 3:numViews
% %     Ps(i,:,:) = algebraicResection(X_init, images{i});
% % end


% Add noise to X_init
X_init = X_init + 0.3*[randn(size(X_init,1)-1, size(X_init,2)); zeros(1,size(X_init,2))];

% Add noise to Ps
% Ps = Ps + 0.001*randn(size(Ps));


%% Bundle Adjustment

% Compute initial reprojection error
reprojErr = computeReprojectionError(Ps, X_init, images);
fprintf('Initial reprojection error: %f\n', norm(reprojErr));

% Perform projective Bundle Adjustment (solves for the projection matrices
% and 3D points)
Ps_opt = Ps;
X_opt = X_init;
[Ps_opt, X_opt, err_opt] = perspectiveBundleAdjustment(Ps_opt, X_opt, images);

fprintf('Reprojection error after Bundle Adjustment: %f\n', err_opt);

% Plot the initial and final cube and camera trajectories
figure;
scatter3(X_init(1,:), X_init(2,:), X_init(3,:), 'filled');
title('Initial structure and camera motion');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis('equal');

figure;
scatter3(X_opt(1,:), X_opt(2,:), X_opt(3,:), 'filled');
title('Optimized structure and camera motion');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis('equal');
