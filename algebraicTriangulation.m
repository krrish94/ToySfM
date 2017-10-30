function X = algebraicTriangulation(x1, x2, P1, P2)
% ALGEBRAICTRIANGULATION  Takes in two sets of correspondences x1, x2, each
% of size 3-by-N (homogeneous image coordinates). Also takes in the
% corresponding 3 x 4 projection matrices P1 and P2. Outputs the
% triangulated 3D points. Does 'vanilla' algebraic triangulation. The
% output points are in 3D homogeneous coordinates, i.e., the output matrix
% is 4 x N.


% Check if number of matches is consistent
if size(x1, 2) ~= size(x2, 2)
    error('x1 and x2 must have the same number of image points.');
end

% Check if the coordinates are homogeneous image coordinates
if size(x1,1) ~= 3 || size(x2,1) ~= 3
    error('x1 and x2 must be homogeneous image coordinates');
end

% % Preprocession (not sure why?) (but, after this results are identical to
% % those from torr_triangulate)
% x1 = x1 ./ 3;
% x2 = x2 ./ 3;

% Initialize the output matrix
X = zeros(4,size(x1,2));

% For each point to be triangulated
for i = 1:size(x1,2)
    
    % Form the linear system JX = 0, where J is a 4 x 4 matrix
    % J results from the equations x1 = P1 * X and x2 = P2 * X
    % J has the following four rows (Pi^(jT) is the jth row vector of Pi)
    % J = [P1^(3T)*x1 - P1^(1T)
    %      P1^(3T)*y1 - P1^(2T)
    %      P2^(3T)*x2 - P2^(1T)
    %      P2^(3T)*y2 - P2^(2T)]
    
    J = zeros(4);
    J(1,:) = x1(1,i).*P1(3,:) - P1(1,:);
    J(2,:) = x1(2,i).*P1(3,:) - P1(2,:);
    J(3,:) = x2(1,i).*P2(3,:) - P2(1,:);
    J(4,:) = x2(2,i).*P2(3,:) - P2(2,:);
    
    % Solve the system using least squares (economy SVD)
    [~, ~, V] = svd(J, 0);
    X(:,i) = V(:,4);
    
    % Note that this hasn't been homogenised so that cheirality can be
    % checked for.
    
end


end
