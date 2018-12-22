function [Ps_opt, X_opt, err_opt] = euclideanGdBundleAdjustment(Ps, X, x)
% EUCLIDEANGDBUNDLEADJUSTMENT  Performs Bundle Adjustment using the 4 x 4
% transformation matrices, camera calibration matrix K, 3D points (4 x N), and M images (x) 
% where all points are assumed to be visible in all the images. Estimates the transformation
% matrices and 3D points by minimizing the reprojection error using the gradient descent optimizer.
%
% Note: For real-world, large-scale BA problems gradient descent (without any pre-conditioning) is a
% poor choice and rarely converges.)

% Number of views
numViews = size(Ps,1);

% Number of points
numPoints = size(X,2);

% Compute initial reprojection error
residual = computeReprojectionError(Ps, X, x);
err = norm(residual,2);

% Optimize only if the reprojection error is above a certian tolerance
% level
tolerance = 1e-16;
if err <= tolerance
    fprintf('Reprojection error is already near-zero. BA not needed.\n');
    Ps_opt = Ps;
    X_opt = X;
    err_opt = 0;
    return;
end

% Set up the parameter vector. The parameter vector has 6 variables for 
% each view (se(3)), and 3 variables for each 3D point.

paramVec = zeros(6*numViews + 3*numPoints, 1);
for i = 1:numViews
    [K,R,t] = decomposeCamera(squeeze(Ps(i,:,:)));
    paramVec(6*(i-1)+1:6*i) = computeKsi([R -R*t; 0 0 0 1]);
end
paramVec(6*numViews+1:end) = reshape(X(1:3,:),[],1);

% Evaluate the Jacobian at the current guess
J = computeLieJacobianMotionAndStructure(paramVec, K, numViews, numPoints);

% Gradient
g = J'*residual;

% Learning rate
lr = 1e-2;

% Vector to store errors over time
errStore = [norm(residual,2)];

% Maximum number of GD iterations
maxIters = 500;

% Stopping criterion (if the gradient is already too small, stop)
stop = (norm(g,'inf') < tolerance);

% Gradient-Descent Iterations
for k = 1:maxIters
    
    if ~stop
        % Update the parameters
        paramVec = paramVec - (lr*g);
        [Ps_new, X_new] = unvectorizeLieParameters(paramVec, K, numViews, numPoints);
            
        % Error resulting from this update
        residual = computeReprojectionError(Ps_new, X_new, x);
        err = norm(residual,2);
        % Store this error for debugging
        errStore = [errStore, err];

        % Update the Jacobian, gradient
        J = computeLieJacobianMotionAndStructure(paramVec, K, numViews, numPoints);
        g = J'*residual;
            
        % Determine if stopping criteria is attained
        if(norm(g,'inf') <= tolerance) || (err <= tolerance)
        	stop = true;
        end
    else
        fprintf('Stopping criteria reached.\n');
        k = k-1;
        break;
    end
    
end

fprintf('Terminated Gradient-Descent iterations. %d iterations complete.\n', k);

% Prepare outputs to be returned by this function
[Ps_opt, X_opt] = unvectorizeLieParameters(paramVec, K, numViews, numPoints);
err_opt = err;

end
