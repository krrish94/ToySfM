function [R,t,matches,P1,P2]=compute_relative_transformation(im1,im2,K)
%
%
%   Compute the relative transformation between the two images. 
%   
%   Input: 
%        im1: image co-ordinates for the first image
%        im2: image co-ordinates for the second image
%          K: intrinsic calibration matrix
%
%   Output:
%          R: Relative rotation between the two images
%          t: Relative translation between the two images
%    matches: The corresponding matches
%         P1: First projection matrix
%         P2: Second projection matrix
%   

    homo_im1=[im1;ones(1,size(im1,2));];    %   homogenize the image co-ordinates
    homo_im2=[im2;ones(1,size(im2,2));];    %   
    
    F=fundmatrix(homo_im1,homo_im2);        %   compute the fundamental matrix
    matches=[homo_im1(1,:)' homo_im1(2,:)' homo_im2(1,:)' homo_im2(2,:)'];
    f=torr_estimateF( matches,1,[50,10],'mapsac',true);
    F=[f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9);];
    F=F';
    E=K'*F*K;                               %   compute the essential matrix
    [u,s,v]=svd(E);
    s(2,2)=s(1,1);   s(3,3)=0;              %   Ensure the rank-2 condition is maintained for the intrinsic calibration matrix
    E=u*s*v';                               %   Reforming the rank-2 essential matrix
    
    matches=[im1(1,:)' im1(2,:)' im2(1,:)' im2(2,:)'];  %   Computing the image matches
    
    [P1,P2,R,t]=torr_linear_EtoPX(E,matches,K,3);       %   Computing the relative transformation via camera resection.
end