function uvPlaneIDData= sr_draw_plane_id(planeProbAcc)

% SC_DRAW_PLANE_ID
%
% Input: 
%   - planeProbAccData
% Output:
%   - uvPlaneIDData

% numPlane = size(planeProbAcc, 1) - 1;
% 
% randSample = rand(1, numUvPix);
% uvPlaneIDData = zeros(1, numUvPix, 'uint8');
% 
% for indPlane = 1: numPlane
%     indSamplePlane = (planeProbAcc(indPlane,:) < randSample ) & ...
%         (planeProbAcc(indPlane + 1, :) >= randSample);
%     uvPlaneIDData(indSamplePlane) = indPlane;
% end

numUvPix = size(planeProbAcc, 1);
numPlane = size(planeProbAcc, 2) - 1;

randSample    = rand(numUvPix, 1);
uvPlaneIDData = zeros(numUvPix, 1, 'uint8');

for indPlane = 1: numPlane
    indSamplePlane = (planeProbAcc(:,indPlane) < randSample ) & ...
        (planeProbAcc(:, indPlane + 1) >= randSample);
    uvPlaneIDData(indSamplePlane) = indPlane;
end


end