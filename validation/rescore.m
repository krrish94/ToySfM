% RESCORE - Repeat re-estimation and re-scoring. To be used by
% RESAMPLE. The process is repeated until either the set does not
% change anymore or maximum maxit time as the set might wobble.
%
% Usage:
% 
% [M, inl, ntinliers]= rescore(x, t, inl, fittingfn, distfn, maxit, low)
%
% Arguments:
%     x         - Data sets to which we are seeking to fit a model M
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
%
%     t         - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%                 Is used in all sampling, resampling and rescoring.
%
%     inl       - Inliers to re-estimate and re-score.
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
%     maxDataTrials - Maximum number of attempts to select a non-degenerate
%                     data set. This parameter is optional and defaults to 100.
%
%     maxTrials - Maximum number of iterations. This parameter is optional and
%                 defaults to 1000.
%     maxit     - Do not repeat rescoring more than maxit times. 
%                 Rescoring is performed until the set does not change
%                 anymore but if the set wobbles maxit prevents it from
%                 getting stuck.
%     low       - How many inliers are required before optimization is done.
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inl       - An array of indices of the elements of x that were
%                 the inliers for the best model.
%     ntinliers - Number of inliers = length(inl)
%
% Copyright (c) 2011-2013 Anders Hast
% Uppsala University
% http://www.cb.uu.se/~aht
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
% AHT 14/1 2013. Extracted from optimalRansac as a function of its own.
% 
%


function [M, inl, ntinliers]= rescore(x, t, inl, fittingfn, distfn, maxit, low)
    j=0;
    npinliers = length(inl);
    pinl=inl;
   
    while j<maxit
        j=j+1;
        M = feval(fittingfn, x(:,inl));
        inl = feval(distfn, M, x, t);
        ntinliers = length(inl);
      
        if ntinliers>low 
            % Is the set changing?
            if ntinliers ~= npinliers
                npinliers=ntinliers;
                pinl=inl;
            elseif sum(pinl==inl)==ntinliers
                % The set is not changing after
                % re-estimation. We are done!
                j=maxit;
            else
                npinliers=ntinliers;
                pinl=inl;
            end
        else
            j=maxit;
        end
    end
end