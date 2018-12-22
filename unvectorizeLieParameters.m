function [Ps, X] = unvectorizeLieParameters(paramVec, K, M, N)
% UNVECTORIZEPARAMETERS  Takes in a parameter vector ((6*M + 3*N) x 1) and
% 'unvectorizes' it, i.e., converts the vector to a set of projection
% matrices and 3D points. K is the camera calibration matrix. M is the number of views,
% and N is the number of 3D points. Note that we return 3D points in homogeneous coordinates,
% i.e., as 4 x 1 vectors, to be consistent with the initialization of 3D
% points, which is also in homogeneous coordinates.

% Projection matrices
Ps = zeros(M,3,4);

for i = 1:M
    tempVar = computeT(paramVec(6*(i-1)+1:6*i));
    Ps(i,:,:) = K * tempVar(1:3,:);
end

% 3D points
X = [reshape(paramVec(6*M+1:end), 3, N); ones(1,N)];

end