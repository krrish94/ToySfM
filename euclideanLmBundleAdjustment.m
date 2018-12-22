function [Ps_opt, X_opt, err_opt] = euclideanLmBundleAdjustment(Ps, X, x)
% EUCLIDEANLMBUNDLEADJUSTMENT  Performs Bundle Adjustment using the 4 x 4
% transformation matrices, camera calibration matrix K, 3D points (4 x N), and M images (x) 
% where all points are assumed to be visible in all the images. Estimates the transformation
% matrices and 3D points by minimizing the reprojection error using the LM optimizer.

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

% Coefficient matrix
A = J'*J;

% Gradient
g = -J'*residual;

% Initialize the damping parameter
dampingCoeff = 10e-2 * max(diag(J'*J));

% Vector to store errors over time
errStore = [norm(residual,2)];

% Number of successful and unsuccessful steps
numSuccess = 0;
numFail = 0;

% Maximum number of LM iterations
maxIters = 500;

% Stopping criterion (if the gradient is already too small, stop)
stop = (norm(g,'inf') < tolerance);

% Levenberg-Marquardt Iterations
for k = 1:maxIters
    
    % Iterate until a feasible point is found
    if ~stop
        % Solve the normal equations for the LM linear system
        deltap = (A + dampingCoeff*eye(size(A,1))) \ g;
        
        % Check if the parameter update is less than tolerance
        if norm(deltap,2) < tolerance * norm(paramVec,2)
            stop = true;
        % If it is not, then compute the updated vector (do not update in
        % the original vector, as we will first determine whether or not it
        % decreases the cost).
        else
            % Updated parameter vector
            paramVec_new = paramVec + deltap;
            
            % Update the parameters (not in the parameter vector, but as
            % projection matrices and 3D points). This function is for
            % easier interpretability, and for computing reprojection
            % error.
            [Ps_new, X_new] = unvectorizeLieParameters(paramVec_new, K, numViews, numPoints);
            
            % Error resulting from this update
            residual_new = computeReprojectionError(Ps_new, X_new, x);
            err_new = norm(residual_new,2);
            
            % Check if the new error is less than the previous error by a
            % margin greater than the tolerance
            if err_new < err
                % Store this error
                errStore = [errStore, err_new];
                if abs(err_new - err) < tolerance
                    stop = true;
                else
                    % This step is successful
                    numSuccess = numSuccess + 1;
                    % Update the parameter vector, Jacobian, and the
                    % residuals, etc.
                    paramVec = paramVec_new;
                    J = computeLieJacobianMotionAndStructure(paramVec, K, numViews, numPoints);
                    residual = residual_new;
                    err = norm(residual,2);
                    A = J'*J;
                    g = -J'*residual;
                    % Determine if stopping criteria is attained
                    stop = (norm(g,'inf') <= tolerance) || (err <= tolerance);
                    % Decrease damping coefficient (since we are in the
                    % trust-region)
                    dampingCoeff = dampingCoeff / 2;
                end
            else
                % We are not in a trust-region
                numFail = numFail + 1;
                % Increase the damping coefficient
                dampingCoeff = dampingCoeff * 2;
            end
        end
    else
        fprintf('Stopping criteria reached.\n');
        k = k-1;
        break;
    end
    
end

fprintf('Terminated Levenberg-Marquardt iterations. %d iterations complete.\n', k);

% Prepare outputs to be returned by this function
[Ps_opt, X_opt] = unvectorizeLieParameters(paramVec, K, numViews, numPoints);
err_opt = err;

end
