function [J] = computeJacobianMotionAndStructure(paramVec, M, N)
% COMPUTEJACOBIANSTRUCTUREANDMOTION  Computes the Jacobian of the
% reprojection error with respect to the motion parameters (projection
% matrices) and the structure parameters (3D points). If there are M views
% and N points, the Jacobian matrix is of size 2 * M * N x (12*M + 3*N).
% It can be viewed as two blocks (left and right blocks), where the left
% block (of size 2 * M * N x 12*M) corresponds to the Jacobian of the
% residual with respect to the motion parameters (projection matrices),
% and the right block (of size 2 * M * N x 3*N) corresponds to the
% Jacobian of the residual with respect to the structure parameters (3D
% points).


% Jacobian matrix initialization
J = zeros(12*M + 3*N, 12*M + 3*N);

% Compute Jacobian with respect to motion and structure parameters
for i = 1:M
    % Current projection matrix guess
    curP = paramVec((i-1)*12+1:12*i);
    % Per view, there will be a 2 x 12 sub-jacobian for each 3D point
    for j = 1:N
        % Current 3D point guess
        curX = [paramVec(12*M + 3*(j-1)+1 : 12*M + 3*j); 1];
        % Defining temporary variables to efficiently compute the
        % sub-jacobian.
        P1X = curP([1,4,7,10])'*curX;
        P2X = curP([2,5,8,11])'*curX;
        P3X = curP([3,6,9,12])'*curX;
        
        % Jacobian with respect to motion parameters
        
        % Computing the sub-jacobian (2 x 12 matrix)
        subjacobian_row1 = [curX' ./ P3X, zeros(1,4), -(P1X/(P3X^2))*curX'];
        subjacobian_row2 = [zeros(1,4), curX' ./ P3X, -(P2X/(P3X^2))*curX'];
        % Insert the sub-jacobian block into its proper position in the
        % Jacobian matrix
        J(2*N*(i-1)+2*(j-1)+1 : 2*N*(i-1)+2*j, 12*(i-1)+1 : 12*i) = ...
            [subjacobian_row1; subjacobian_row2];
        clear subjacobian_row1 subjacobian_row2
        
        % Jacobian with respect to structure parameters
        
        % Computing the sub-jacobian (2 x 3 matrix)
        subjacobian_row1 = (1/P3X)^2 * [ curP(1)*P3X - P1X*curP(3), ...
            curP(4)*P3X - P1X*curP(6), curP(7)*P3X - P1X*curP(9)];
        subjacobian_row2 = (1/P3X)^2 * [curP(2)*P3X - P2X*curP(3), ...
            curP(5)*P3X - P2X*curP(6), curP(8)*P3X - P2X*curP(9)];
        % Insert the sub-jacobian block into its proper position in the
        % Jacobian matrix
        J(2*N*(i-1)+2*(j-1)+1 : 2*N*(i-1)+2*j, 12*M+3*(j-1)+1 : 12*M+3*j) = ...
            [subjacobian_row1; subjacobian_row2];
        clear subjacobian_row1 subjacobian_row2
        
    end
end

end
