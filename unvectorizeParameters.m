function [Ps, X] = unvectorizeParameters(paramVec, M, N)
% UNVECTORIZEPARAMETERS  Takes in a parameter vector ((12*M + 3*N) x 1) and
% 'unvectorizes' it, i.e., converts the vector to a set of projection
% matrices and 3D points. M is the number of views, and N is the number of
% 3D points. Note that we return 3D points in homogeneous coordinates,
% i.e., as 4 x 1 vectors, to be consistent with the initialization of 3D
% points, which is also in homogeneous coordinates.


% Projection matrices
Ps = zeros(M,3,4);

for i = 1:M
    Ps(i,:,:) = reshape(paramVec(12*(i-1)+1:12*i), 3, 4);
end

% 3D points
X = [reshape(paramVec(12*i+1:end), 3, N); ones(1,N)];

end
