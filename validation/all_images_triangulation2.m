function [reconstructed_3d_pts] = all_images_triangulation2(Rs,ts,K,images)
%
%   Input:  
%       Rs: Rotations in a single co-ordinate frame
%       ts: translations in a single co-ordinate frame
%   images: cell of image co-ordinates
%   
%   Output: 
%   pts_3d_reconstructed: Reconstructed 3d points.    
%   


num_pts=size(images{1,1},2);    %   Number of image points
num_cameras=size(Rs,1);          %   Number of cameras       %   


point_tracks=zeros(2*num_pts,num_cameras);  %   Re-organizing the images

for i=1:num_pts
    for j=1:num_cameras
        point_tracks(2*(i-1)+1:2*(i-1)+2,j)=images{j}(1:2,i);
    end
end

reconstructed_3d_pts=zeros(3,num_pts);

for i=1:num_pts
    track=point_tracks(2*(i-1)+1:2*(i-1)+2,:);
    reconstructed_3d_pts(:,i)=triangulate_algebraic_least_sq(Rs,ts,K,track);
end


end

function [reconstructed_pt]=triangulate_algebraic_least_sq(Rs,ts,K,track)
%
%   Input: 
%       Rs: Rotations
%       ts: Translations
%        K: Intrinsic calibration matrix
%    track: Camera track
%
%   Output: 
%   reconstructed_pt: Reconstructed 3d pt.

num_cameras=size(Rs,1);

A=zeros(2*num_cameras,3);
b=zeros(2*num_cameras,1);

for i=1:num_cameras
    x=track(1,i);   %   x img-co-ordinate
    y=track(2,i);   %   y img-co-ordinate
    R=squeeze(Rs(i,:,:));
    t=ts(:,i);
    kr=K*R';
    kt=K*(-R'*t);
   A(2*(i-1)+1,:)=x*kr(3,:)-kr(1,:);
   A(2*(i-1)+2,:)=y*kr(3,:)-kr(2,:);
   b(2*(i-1)+1,1)=kt(1,1)-x*kt(3,1);
   b(2*(i-1)+2,1)=kt(2,1)-y*kt(3,1);
end

reconstructed_pt=linsolve(A,b);

end