% OPTIMALRANSAC - Robustly fits a model to data with the Optimal RANSAC algorithm
%
% Usage:
%
% [M, sinliers] = optimalRansac(x, fittingfn, distfn, s, t, acc, ...
%                        maxDataTrials, MaxTrials, maxit, nrtrials, low, ner)
%
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
%                 Is used in all sampling, resampling and rescoring.
%     acc       - The same as above but is used to prune the set. If acc<t
%                 then pruning is performed.
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
%     nrtrials  - Defines how many resamplings shall be done.
%     low       - How many inliers are required before optimization is done.
%     ner       - Number of equal sets required before stopping.
%                 1 means that one pair of equal sets have been obtaines.
%                 2 Means that three sets were equal etc.
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
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
% Rewritten from and to be used with code that is Copyright (c) 2003-2006 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au    
% http://www.csse.uwa.edu.au/~pk
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
% AHT 10/6 2015. Send low as parameter to PRUNESET and s as parameter to
% RESAMPLE
% AHT 10/6 2015. Changed computation of iterations
%
function [M, sinliers] = optimalRansac(x, fittingfn, distfn, s, t, acc, ...
                               maxDataTrials, MaxTrials, maxit, nrtrials, low, ner)
                  
    % Test number of parameters
    narginchk ( 7, 12 );
    nargoutchk ( 2, 2 );
    
    if nargin == 7
        maxDataTrials=100;
        MaxTrials=10000;
        maxit=20;
        nrtrials=8;
        low=5;
        ner=1;
    end
                      
    nes=0;      % Number of equal sets
    maxinl=0;   % Maximum number of inliers    
    iter=0;
    
    while nes< ner & iter < MaxTrials
        iter=iter+1;
        [M,inliers, ninliers] = sampleSet(x, fittingfn, distfn, s, t, maxDataTrials);
        
        if ninliers>low       
            [M, inliers, ninliers]= resample(x, M, s, t, inliers, ninliers, fittingfn, distfn, nrtrials, maxit, low);

            %  Clean the set up!
            if acc<t
                % Here we should have something good so prune the set
                % Sometimes the set might have one inlier too much, perhaps
                % because of rounding errors. The tactic below seem to
                % efficiently avoid such situations.
                [M, putative, ninliers]= pruneset(M, x(:,inliers), acc, low, fittingfn, distfn);
                inliers=inliers(putative);
            end
        end

        % Did we get the same set twice?
        if maxinl>low && maxinl==ninliers
            if sum(inliers==sinliers)==maxinl % Yes we got the same set twice
                nes=nes+1;
            else % Not the same set 
                nes=0;
                sinliers=inliers;
            end
        elseif ninliers>maxinl 
            nes=0;
            maxinl=ninliers;
            sinliers=inliers;
        elseif ninliers==maxinl-1
            nes=0;
            maxinl=ninliers;
            sinliers=inliers;
        end       
        
    end
end
