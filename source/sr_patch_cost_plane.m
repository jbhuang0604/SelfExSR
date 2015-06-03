
function costPlane = sr_patch_cost_plane(mLogLPlaneProb, uvPlaneIDData, trgPixInd, srcPixSub)

% SR_PATCH_COST_PLANE
% 
% Compute planar costs (See Eqn 11 in the paper)
% 
% Jia-Bin Huang, Sing Bing Kang, Narendra Ahuja, Johannes Kopf
% Image Completion using Planar Structure Guidance
% ACM Transactions on Graphics (Proceedings of SIGGRAPH 2014).
% 
% Input:
%   - mLogLPlaneProb: minus log-likelihood of plane assignment
%   - uvPlaneIDData:  plane ID for each target patch
%   - trgPixInd:      target patch index
%   - srcPixSub:      source patch position
% Output:
%   - costPlane:      plane compatibility cost
% =========================================================================

[imgH, imgW, numPlane] = size(mLogLPlaneProb);

% Get the source patch index
srcPixSub = round(srcPixSub);

srcPixSub(:,1) = sr_clamp(srcPixSub(:,1), 1, imgW);
srcPixSub(:,2) = sr_clamp(srcPixSub(:,2), 1, imgH);

srcPixInd = sub2ind([imgH, imgW, numPlane], srcPixSub(:,2), srcPixSub(:,1), single(uvPlaneIDData));

% Plane compatibility cost
costPlane = mLogLPlaneProb(trgPixInd) + mLogLPlaneProb(srcPixInd);

end 