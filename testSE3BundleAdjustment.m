%% Script to test the ToySfM example

% The scene and camera trajectories are similar to those in the simple
% example script in 'testProjectiveBundleAdjustment.m'. For comments, one
% may refer to that script. This script has very few comments for the scene
% construction part.


%% Synthesize the scene

% Seed the random number generator (for repeatability)
rng('default');

% Camera intrinsics
K = eye(3);
% Ground-Truth world points
worldPoints_gt = generateCube(4,16);
% Number of views
numViews = 8;
% Number of points
numPoints = size(worldPoints_gt,2);
% Simulate the camera trajectory (circle) around the cube
[Rs, ts] = generateCameraTrajectory(numViews, 25);
% Camera locations
camLocs = ts';
% Initialize a variable to hold all projection matrices
Ps = zeros(numViews, 3, 4);

% % Visualize the scene and the camera trajectories
% scatter3(worldPoints_gt(1,:), worldPoints_gt(2,:), worldPoints_gt(3,:), 'filled');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('The synthetic scene and cameras');
% axis('equal');
% hold on;
% for i = 1:numViews
%     R_cam_W = squeeze(Rs(i,:,:));
%     % Plot the camera frustum
%     plotCameraFrustum(squeeze(Rs(i,:,:)), camLocs(:,i));
%     % Plot the camera coordinate frame
%     scatter3(camLocs(1,i), camLocs(2,i), camLocs(3,i), 25, [0, 0, 1], 'filled');
%     text(camLocs(1,i), camLocs(2,i), camLocs(3,i)-2, sprintf('C%d',i), 'Color', 'red');
%     quiver3(camLocs(1,i), camLocs(2,i), camLocs(3,i), R_cam_W(1,1), R_cam_W(2,1), R_cam_W(3,1), 'r', 'LineWidth', 2);
%     text(camLocs(1,i) + R_cam_W(1,1), camLocs(2,i) + R_cam_W(2,1), camLocs(3,i) + R_cam_W(3,1), 'X', 'Color', 'red');
%     quiver3(camLocs(1,i), camLocs(2,i), camLocs(3,i), R_cam_W(1,2), R_cam_W(2,2), R_cam_W(3,2), 'g', 'LineWidth', 2);
%     text(camLocs(1,i) + R_cam_W(1,2), camLocs(2,i) + R_cam_W(2,2), camLocs(3,i) + R_cam_W(3,2), 'Y', 'Color', 'green');
%     quiver3(camLocs(1,i), camLocs(2,i), camLocs(3,i), R_cam_W(1,3), R_cam_W(2,3), R_cam_W(3,3), 'b', 'LineWidth', 2);
%     text(camLocs(1,i) + R_cam_W(1,3), camLocs(2,i) + R_cam_W(2,3), camLocs(3,i) + R_cam_W(3,3), 'Z', 'Color', 'blue');
%     clear R_cam_W;
% end



% Synthesize images of the cube
images = synthesizeImages(worldPoints_gt, K, Rs, ts);



%% Construct the initial guess


% Compute relative poses of each view with respect to the first view
for i = 2:numViews
    
    [R_21, t_21] = computeRelativePose(images{1}, images{i}, K);
    
    % Initialize structure
    if i == 2
        
        Ps(1,:,:) = K * [eye(3), [0; 0; 0]];
        Ps(2,:,:) = K * [R_21, t_21];
        X_init = algebraicTriangulation(images{1}, images{2}, squeeze(Ps(1,:,:)), squeeze(Ps(2,:,:)));
        X_init = X_init ./ repmat(X_init(4,:), 4, 1);
       
    % For all other views
    else
        
        % Compute the scaling factor and scale the translation estimate
        P_temp = K * [R_21, t_21];
        X_temp = algebraicTriangulation(images{1}, images{i}, squeeze(Ps(1,:,:)), P_temp);
        X_temp = X_temp ./ repmat(X_temp(4,:), 4, 1);
        
        numSamplesForScale = 25;
        scaleFactors = zeros(numSamplesForScale,1);
        for j = 1:numSamplesForScale
            sampledPoints = randsample(numPoints,2);
            d1 = sqrt(sum((X_init(1:3,sampledPoints(1)) - X_init(1:3,sampledPoints(2))).^2));
            di = sqrt(sum((X_temp(1:3,sampledPoints(1)) - X_temp(1:3,sampledPoints(2))).^2));
            scaleFactor = d1 / di;
            scaleFactors(j) = scaleFactor;
        end
        avgScaleFactor = mean(scaleFactors);
        
        t_21 = avgScaleFactor * t_21;
        
        % Construct the projection matrix using the estimated relative pose
        Ps(i,:,:) = K * [R_21, t_21];
        
    end
    
    % Clear temporary variables
    clear R_21 t_21
    
end


% Add noise to X_init
X_init = X_init + 0.3*[randn(size(X_init,1)-1, size(X_init,2)); zeros(1,size(X_init,2))];

% % Add noise to Ps
% Ps = Ps + 0.01*randn(size(Ps));


%% Bundle Adjustment on the SE(3) manifold
% We do it over the Lie algebra se(3), actually!

% Compute initial reprojection error
reprojErr = computeReprojectionError(Ps, X_init, images);
fprintf('Initial reprojection error: %f\n', norm(reprojErr));

% Perform projective Bundle Adjustment (solves for the projection matrices
% and 3D points)
Ps_opt = Ps;
X_opt = X_init;
[Ps_opt, X_opt, err_opt] = perspectiveBundleAdjustment(Ps_opt, X_opt, images);

fprintf('Reprojection error after Bundle Adjustment: %f\n', err_opt);

% % Plot the initial and final cube and camera trajectories
% figure;
% scatter3(X_init(1,:), X_init(2,:), X_init(3,:), 'filled');
% title('Initial structure and camera motion');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis('equal');
% 
% figure;
% scatter3(X_opt(1,:), X_opt(2,:), X_opt(3,:), 'filled');
% title('Optimized structure and camera motion');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis('equal');