function uvValidInd = sr_check_valid_uv(srcPos, validSrcMask)

uvSub = round(srcPos);

% clamp
uvSub(1,:) = sr_clamp(uvSub(1,:), 1, size(validSrcMask,2));
uvSub(2,:) = sr_clamp(uvSub(2,:), 1, size(validSrcMask,1));

% uvSub = round(uvSub);

uvInd = sub2ind(size(validSrcMask), uvSub(2,:), uvSub(1,:));

uvValidInd = validSrcMask(uvInd);

% uvValid.ind(uvValid.ind) = uvValid.ind(uvValid.ind) & uvValidCand;
% uvValid.pos = find(uvValid.ind);

end