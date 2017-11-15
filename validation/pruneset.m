% PRUNESET -  Try to obtain the best homography from the whole set of points
% using the already computed homography, which can be computed by any
% RANSAC alhorithm. It is called by OPTIMALRANSAC.
%
% Usage:
%
% [H, putative, len]= pruneset(H, x, acc, fittingfn, distfn)
%
%
% Arguments:
%     H         - The model Homography for x.
%
%     x         - Data sets to which we are seeking to fit a model H
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
%     acc       - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%                 Is used in all sampling, resampling and rescoring.
%
%     fittingfn - Handle to a function that fits a model to s
%                 data from x.  It is assumed that the function is of the
%                 form: 
%                    M = fittingfn(x)
%                 Note it is possible that the fitting function can return
%                 multiple models (for example up to 3 fundamental matrices
%                 can be fitted to 7 matched points).  In this case it is
%                 assumed that the fitting function returns a cell array of
%                 models.
%                 If this function cannot fit a model it should return M as
%                 an empty matrix.
%
%     distfn    - Handle to a function that evaluates the
%                 distances from the model to data x.
%                 It is assumed that the function is of the form:
%                    [inliers, M] = distfn(M, x, t)
%                 This function must evaluate the distances between points
%                 and the model returning the indices of elements in x that
%                 are inliers, that is, the points that are within distance
%                 't' of the model.  Additionally, if M is a cell array of
%                 possible models 'distfn' will return the model that has the
%                 most inliers.  If there is only one model this function
%                 must still copy the model to the output.  After this call M
%                 will be a non-cell object representing only one model. 
%
% Returns:
%     H         - The model for the pruned set.
%     putative  - An array of indices of the elements of x that were
%                 the inliers for the pruned set.
%     iterations - How many iterations was done. This is measuered in 
%                 terms of ordinary RANSAC iterations so it can be 
%                 compared to RANSAC. 
%
% For an example of the use of this function see OPTIMALRANSAC 
%
% Copyright (c) 2011-2013 Anders Hast
% Uppsala University
% http://www.cb.uu.se/~aht
% 
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% History
% AHT 30/1 2013. homogdist2d in OPTIMALRANSACFITHOMOGRAPHY was changed to 
%                return also the distances. This is not done in the original
%                code by Peter Kovesi. 
% AHT 10/6 2015. Use low instead of hardcoded value. Also check that FEVAL
% will only be called if not too many are removed.

function [H, putative, len]= pruneset(H, x, acc, low, fittingfn, distfn)
    len=size(x,2);
    
    repeat=1;
    putative=[1:len];
    while len>low & repeat
        % Compute the distances
        [inliers, d2] = feval(distfn, H, x(:,putative), acc);
        
        % Remove the most extreme point
        [t1,rc]=max(d2);
        if t1>acc & len-length(rc)> low
            putative(rc)=[];       
            H = feval(fittingfn, x(:,putative));
        else
            repeat=0;
        end
  
        len=size(putative,2);
    end
end