function map = sr_update_uvMap(map, data, uvPixInd)

% SC_UPDATE_UVMAP
%
% Update the map from the data and the indices

[imgH, imgW, nCh] = size(map);

offset = uint64((0:nCh-1)*imgH*imgW);

uvPixInd = bsxfun(@plus, uvPixInd, offset);
map(uvPixInd) = data;

end