% SAMPLESET - Sample a non degenerate set. To be used by OPTIMALRANSAC.
% This code was extraced with minor changes from RANSAC.
%
% Usage:
%
% [M, inliers, ninliers] = sampleSet(x, fittingfn, distfn, s, t, maxDataTrials)
%
% Arguments:
%     x         - Data sets to which we are seeking to fit a model M
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
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
%     s         - The minimum number of samples from x required by
%                 fittingfn to fit a model.
%
%     t         - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%
%     maxTrials - Maximum number of iterations. This parameter is optional and
%                 defaults to 1000.
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
%     ninliers  - Number of inliers.
%
%
% For an example of the use of this function see RANSACFITHOMOGRAPHY or
% RANSACFITPLANE 
%
% References:
%    M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
%    for model fitting with applications to image analysis and automated
%    cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981
%
%    Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
%    Computer Vision". pp 101-113. Cambridge University Press, 2001
%
% The original code:
%   Copyright (c) 2003-2006 Peter Kovesi
% Changes:
%   The code was extracted from RANSAC by Peter Kovesi and some small
%   cchanges were made.
%   Copyright (c) 2011-2013 Anders Hast
%
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au    
% http://www.csse.uwa.edu.au/~pk
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

function [M,inliers, ninliers] = sampleSet(x, fittingfn, distfn, s, t, maxDataTrials)
                  
    % Test number of parameters
    narginchk ( 6,6 );
    nargoutchk ( 3, 3 );
        
    [rows, npts] = size(x);  
        
    % Select at random s datapoints to form a trial model, M.
    % In selecting these points we have to check that they are not in
    % a degenerate configuration.
    degenerate = 1;
    count = 0;  
    
    while degenerate
        % Generate s random indicies in the range 1..npts
        % (If you do not have the statistics toolbox, or are using Octave,
        % use the function RANDOMSAMPLE from my webpage)
        % ind = randsample(npts, s);

        ind = randomsample(npts, s);
        % Assume we have a degenerate case
        degenerate=0;
 
        % Fit model to this random selection of data points.
        % Note that M may represent a set of models that fit the data in
        % this case M will be a cell array of models
        M = feval(fittingfn, x(:,ind));
      
        % Depending on your problem it might be that the only way you
        % can determine whether a data set is degenerate or not is to
        % try to fit a model and see if it succeeds.  If it fails we
        % reset degenerate to true.
        if sum(sum(M==0))==9 %isempty(M)
            degenerate = 1;
        end
           
        % Safeguard against being stuck in this loop forever
        count = count + 1;
        if count > maxDataTrials
            warning('Unable to select a nondegenerate data set');
            break
        end
    end
        
    % Once we are out here we should have some kind of model...        
    % Evaluate distances between points and model returning the indices
    % of elements in x that are inliers.  Additionally, if M is a cell
    % array of possible models 'distfn' will return the model that has
    % the most inliers.  After this call M will be a non-cell object
    % representing only one model.
    [inliers] = feval(distfn, M, x, t);

    % Find the number of inliers to this model.
    ninliers = length(inliers); 
end
 