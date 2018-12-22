function [ksi] = computeKsi(T)
% COMPUTEKSI Convert given SE(3) matrix to se(3) vector representation
% using the logarithmic map (Log).
%
% T is a 4 x 4 SE(3) matrix
% ksi is a 6 x 1 se(3) vector

     R = T(1:3,1:3);
     thetha = acos((trace(R) - 1)/2);
     
     % Handling zero rotation
     if(thetha == 0)
         w = [0;0;0];
         v = T(1:3,4);
         ksi = [v; w];
     else
         w = thetha * (1/(2*sin(thetha))*[R(3,2) - R(2,3);R(1,3) - R(3,1);R(2,1) - R(1,2)]);
         wx = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0;];
         t = T(1:3,4);
         % Computing the linear velocity 3 x 1 vector (expression from ethan eade doc)
         v = (eye(3) - (1/2 * (wx)) + (((1/(thetha * thetha)) * (1 - ((thetha * sin(thetha)) / (2*(1 - cos(thetha)))))) * (wx * wx)))*t;
         ksi = [v; w];
     end     
end
