%% Temporary script used to test parts of code


%% Test 1: Coordinate frame conventions for Essential matrix decomposition

% % Check the convention of the R's t's with those estimated from the
% % fundamental matrix
% 
% R1 = squeeze(Rs(1,:,:));
% t1 = ts(1,:)';
% R2 = squeeze(Rs(2,:,:));
% t2 = ts(2,:)';
% % Transform that maps points from camera coordinates of first view to the
% % camera coordinates of the second view
% mat1 = inv([R2, t2; zeros(1,3), 1]) * [R1, t1; zeros(1,3), 1];
% % Obtained from fundamental matrix decomposition
% mat2 = [R_21, t_21; zeros(1,3), 1];
% 
% % We now verify that mat1 and mat2 perform the same task, i.e., transfer
% % points from the first camera frame to the second camera frame. We verify
% % that by projection 3D structure estimated from triangulating the first
% % two views onto the second image.
% 
% % Transform matrix (4 x 4) that transforms from cam 1 coordinates to cam 2
% % coordinates
% pts_cam2 = mat2 * X_init;
% pts_cam2 = pts_cam2 ./ repmat(pts_cam2(4,:),4,1);
% pts_img2 = K * pts_cam2(1:3,:);
% pts_img2 = pts_img2 ./ repmat(pts_img2(3,:),3,1);
% % Compute reprojection error
% err_tmp = abs(pts_img2 - images{2});
% err_tmp(find(err_tmp < 1e-5)) = 0;
% err_tmp = sum(sum(err_tmp));
% clear pts_cam2 pts_img2
% 
% % For mat1, we need to use the 'real-sized' cube, expressed in cam1
% % coordinates.
% pts_cam1 = inv([R1 t1; zeros(1,3), 1]) * [worldPoints_gt; ones(1,size(worldPoints_gt,2))];
% pts_cam2 = mat1 * pts_cam1;
% clear pts_cam1
% pts_cam2 = pts_cam2 ./ repmat(pts_cam2(4,:),4,1);
% pts_img2 = K * pts_cam2(1:3,:);
% pts_img2 = pts_img2 ./ repmat(pts_img2(3,:),3,1);
% % Compute reprojection error
% err_tmp = abs(pts_img2 - images{2});
% err_tmp(find(err_tmp < 1e-5)) = 0;
% err_tmp = sum(sum(err_tmp));
% clear pts_cam2 pts_img2

% clear mat1 mat2 R1 R2 t1 t2


%% Test 2: Testing View 3 onwards

% % Transforms that take a point from the camera coordinates to world
% % coordinates (for camera i, the Transform is Ti)
% T1 = [squeeze(Rs(1,:,:)), ts(1,:)'; zeros(1,3), 1];
% T2 = [squeeze(Rs(2,:,:)), ts(2,:)'; zeros(1,3), 1];
% T3 = [squeeze(Rs(3,:,:)), ts(3,:)'; zeros(1,3), 1];
% T4 = [squeeze(Rs(4,:,:)), ts(4,:)'; zeros(1,3), 1];
% T5 = [squeeze(Rs(5,:,:)), ts(5,:)'; zeros(1,3), 1];
% 
% % Transforms from camera frame 1 to other camera coordinates
% T_12 = inv(T1) * T2;
% T_13 = inv(T1) * T3;
% % T_14 = inv(T1) * T4;
% % T_15 = inv(T1) * T5;
% 
% clear T1 T2 T3 T4 T5
% 
% % Transforms from world coordinates to camera i (estimated via fundamental
% % matrix decomposition).
% T12_est = [inv(K)*squeeze(Ps(2,:,:)); zeros(1,3), 1];
% T13_est = [inv(K)*squeeze(Ps(3,:,:)); zeros(1,3), 1];
% % T14_est = [inv(K)*squeeze(Ps(4,:,:)); zeros(1,3), 1];
% % T15_est = [inv(K)*squeeze(Ps(5,:,:)); zeros(1,3), 1];
% 
% % Verify that these sets of matrices are inverses of each other, i.e., that
% % they compose to the identity element of SE(3).
% T_12*T12_est
% T_13*T13_est
% % T_14*T14_est
% % T_15*T15_est


%% Test 3: Test image synthesis

% % View to test
% i = 2;
% 
% % Camera coordinates
% R = squeeze(Rs(i,:,:));
% t = ts(i,:)';
% X_cam = R' * (worldPoints_gt - t);
% scatter3(X_cam(1,:), X_cam(2,:), X_cam(3,:), 'filled');
% view(0,-90);
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% mean(X_cam,2)
% 
% % Synthesize image
% x_im = K * X_cam;
% x_im = x_im ./ repmat(x_im(3,:), 3, 1);
% figure;
% scatter(x_im(1,:), x_im(2,:), 'filled');

% clear X_cam x_im R t
