function uvTformDataT = sr_src_domain_tform(uvPlaneID, modelPlane, uvTformAData, srcPos, trgPos)
% SR_SRC_DOMAIN_TFORM
%
% The function estimates the source patch domain transformation based on
% the source and target position, the plane parameters and apply the affine
% transformation
%
% Input:
%     - uvPlaneID:    [numUvPix] x [1], the plane label of uv pixels
%     - modelPlane:   plane model
%     - uvTformAData: [numUvPix] x 4, affine transformation
%     - srcPos:       [numUvPix] x [2], the center position of source patch
%     - trgPos:       [numUvPix] x [2], the center position of target patch
% Output:
%     - uvTformDataT: [numUvPix] X [9]
%
% The function of computing homography transformation induced by plane 
% is modified from the following paper:
% 
% Jia-Bin Huang, Sing Bing Kang, Narendra Ahuja, and Johannes Kopf,
% Image Completion using Planar Structure Guidance,
% ACM Transactions on Graphics (Proceedings of SIGGRAPH 2014), 33(4), 2014
% =========================================================================

% TO-DO: Make sure that the predicted transformation can point it 
% the desired source patch position

% =========================================================================
% Compute homography transformation induced by plane
% =========================================================================

numUvPix = size(srcPos, 1);
uvTformDataH = zeros(numUvPix, 9, 'single');
I = eye(3);

for indPlane = 1: modelPlane.numPlane
    % The rectifying transformation for the plane
    rectMat = modelPlane.rectMat{indPlane};
    h7 = rectMat(3,1);
    h8 = rectMat(3,2);
    
    % Retrieve the uv pixels that have the current plane label
    uvPlaneIndCur = uvPlaneID == indPlane;
    numPlanePixCur = sum(uvPlaneIndCur);
    
    if(numPlanePixCur)
        % Target patch center position in the rectified domain
        trgPosCur = trgPos(uvPlaneIndCur, :) - 1;
        trgPosCurR = sr_apply_tform_H(trgPosCur, h7, h8);
        
        % Source patch center position in the rectified domain
        srcPosCur = srcPos(uvPlaneIndCur, :) - 1;
        srcPosCurR = sr_apply_tform_H(srcPosCur, h7, h8);
        
        % Displacement vector from target to source position
        dRect = srcPosCurR - trgPosCurR;
        
        % Compute the transformation that maps from target to source
        % (See Eqn 8 in the paper [Huang et al. TOG 2014])
        uvTformCur = zeros(numPlanePixCur, 9, 'single');
        uvTformCur(:,[1,4,7]) = bsxfun(@times, dRect(:,1), [h7, h8, 1]);
        uvTformCur(:,[2,5,8]) = bsxfun(@times, dRect(:,2), [h7, h8, 1]);
        
        dTemp = dRect*[h7;h8]; % dTemp = dx*h7 + dy*h8
        uvTformCur(:,[3,6,9]) = bsxfun(@times, dTemp, -[h7, h8, 1]);
        uvTformCur = bsxfun(@plus, uvTformCur, I(:)');
        
        % Apply the offset to cancel out the dependency of the target position
        % (See Eqn 9 in the paper [Huang et al. TOG 2014])
        uvTformDataH(uvPlaneIndCur, :)   = sr_trans_tform(uvTformCur, trgPosCur);
        uvTformDataH(uvPlaneIndCur, 7:8) = uvTformDataH(uvPlaneIndCur, 7:8) + 1;
    end
end

% =========================================================================
% Apply the similarity transformation
% =========================================================================
% The transformation T first map from original reference points to affine
% transformation coordinate using A, and then map to the desired target
% position using H
% -> T = H * A
uvTformDataT = uvTformDataH;
uvTformDataT(:,1:3) = bsxfun(@times, uvTformDataH(:,1:3), uvTformAData(:,1)) + ...
    bsxfun(@times, uvTformDataH(:,4:6), uvTformAData(:,2));
uvTformDataT(:,4:6) = bsxfun(@times, uvTformDataH(:,4:6), uvTformAData(:,3)) + ...
    bsxfun(@times, uvTformDataH(:,4:6), uvTformAData(:,4));

uvTformDataT = bsxfun(@rdivide, uvTformDataT, uvTformDataT(:,9));

end


function x = sr_apply_tform_H(x, h7, h8)

% Apply homography H with third row [h7, h8, 1] to 2D points x

y = x(:,1)*h7 + x(:,2)*h8 + 1;
x = bsxfun(@rdivide, x(:,1:2), y + eps);

% numPix = size(x, 1);
% x = cat(2, x, ones(numPix, 1, 'single'));
% x = x*H';
% x = bsxfun(@rdivide, x(:,1:2), x(:,3) + eps);

end