function planeProbAcc = sr_prep_plane_prob_acc(planeProb, uvPixInd)

numPlane = size(planeProb, 3);
numUvPix = size(uvPixInd, 1);

planeProbAcc = zeros(numUvPix, numPlane + 1, 'single');

% Compute the accumulative probability
for i = 1: numPlane
    planeProbCur = planeProb(:,:,i);
    planeProbAcc(:, i+1) = planeProbCur(uvPixInd);
    % Accumulative probability, starting from 0
    if(i ~= 1)
        planeProbAcc(:, i+1) = planeProbAcc(:, i+1) + planeProbAcc(:, i);
    end
end

end