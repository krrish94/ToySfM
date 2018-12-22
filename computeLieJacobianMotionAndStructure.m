function [J] = computeLieJacobianMotionAndStructure(paramVec, K, M, N)
% COMPUTELIEJACOBIANSTRUCTUREANDMOTION  Computes the Jacobian of the
% reprojection error with respect to the motion parameters (se(3) vectors)
% and the structure parameters (3D points).
%
% M - Number of views
% N - Number of features
% Size of paramVec - [(6*M+3*N) x 1]

fx = K(1,1);
fy = K(2,2);

J = zeros(2*N*M, (6*M + 3*N));

for i = 1:M
    
    currXi = paramVec(6*(i-1)+1 : 6*i);
    T = computeT(currXi);
    
    for p = 1:N
        currXp = paramVec(6*M+3*(p-1)+1 : 6*M+3*p);
        currXp = T * [currXp;1];
        x = currXp(1);
        y = currXp(2);
        z = currXp(3);
        
        % Refer C. Kerl's thesis
        Jmotion = [fx * 1/z, 0, -fx * (x/z^2), -fx * (x * y / z*2), fx * (1 + x^2/z^2), -fx * y/z; ...
            0, fy * 1/z, -fy * y / z^2, -fy * (1 + y^2/z^2), fy * x * y / z^2, fy * x/z];
        
        J(2*N*(i-1)+2*(p-1)+1 : 2*N*(i-1)+2*p, 6*(i-1)+1 : 6*i) = Jmotion;
        clear Jmotion
        
        Jstructure = [fx / z, 0, -fx * x / z^2;
            0, fy / z, -fy * y / z^2] * T(1:3,1:3);
        
        J(2*N*(i-1)+2*(p-1)+1 : 2*N*(i-1)+2*p,6*M+3*(p-1)+1 : 6*M+3*p) = Jstructure;
        clear Jstructure
        
    end
end
end