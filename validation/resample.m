% RESAMPLE - Repeat re-sampling and call RESCORE for re-estimation and re-scoring.
% To be used by OPTIMALRANSAC. The process is repeated until either the set does not
% change anymore or maximum maxit time as the set might wobble.
%
% Usage:
% 
% [M, inliers, ninliers]= resample(x, M, t, inliers, ninliers, fittingfn, distfn, nrtrials, maxit, low)
%
% Arguments:
%     x         - Data sets to which we are seeking to fit a model M
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
%
%     M         - The model having the greatest number of inliers.
%
%     t         - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%                 Is used in all sampling, resampling and rescoring.
%
%     inliers   - Inliers to re-sample and re-estimate and re-score.
%
%     ninliers  - Number of inliers.
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
%     nrtrials  - Number if re-samplings. If a better set is found the
%                 algorithm starts all over to do the re-sampliings nrtrials
%                 times on that larger set.
%
%     maxit     - Do not repeat rescoring more than maxit times. 
%                 Rescoring is performed until the set does not change
%                 anymore but if the set wobbles maxit prevents it from
%                 getting stuck.
%     low       - How many inliers are required before optimization is done.
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
%     ninliers - Number of inliers = length(inliers)
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
% AHT 10/6 2015. Use s instead of hardcoded value in computation of set
% size in parameter to RANDOMSAMPLE
%
function [M, inliers, ninliers]= resample(x, M, s, t, inliers, ninliers, fittingfn, distfn, nrtrials, maxit, low)
    % Try using a subset of the set in order to find a
    % better mix to use for matching.        
 
    xx=x(:,inliers);
    i=0;
    while i<nrtrials
        i=i+1;
        ind = randomsample(ninliers, max(s,round(ninliers/4)));
        H = feval(fittingfn, xx(:,ind));
        inl = feval(distfn, H, x, t);

        if length(inl)>low
            [H, inl, ntinliers]= rescore(x, t, inl, fittingfn, distfn, maxit, low);

            % Do we have a better set?!
            if ntinliers>ninliers
                xx=x(:,inl);
                M=H;
                ninliers=ntinliers;
                inliers=inl;
                i=0;
            end
        end
    end
end