function     uvTformScale = sr_scale_tform(H)

% SR_SCALE_TFORM

% Estimate the approximated scale from the homography transformation using
% first-order Tayor expansion.
% 
% For more detail, please see the paper
%
% [Chum and Matas 10] 
% Chum, Ondrej, and Jirí Matas. 
% "Planar affine rectification from change of scale." 
% Asian Conference on Computer Vision, ACCV 2010
%
%      h1 h4 h7       h1 - h7*h3   h4 - h7*h6   h7       1  0  0 
%T = [ h2 h5 h8 ] = [ h2 - h8*h3   h5 - h8*h6   h8 ] * [ 0  1  0 ]
%      h3 h6 h9           0            0         1       h3 h6 1
% T = A*H

% Here the scale estimation is the determinant of the matrix A

uvTformScale = (H(:,1) - H(:,7).*H(:,3)).* (H(:,5) - H(:,8).*H(:,6)) ...
             - (H(:,4) - H(:,7).*H(:,6)).* (H(:,2) - H(:,8).*H(:,3));

uvTformScale =  sqrt(abs(uvTformScale));

end