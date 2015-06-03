function uvTform = sr_trans_tform(uvTformV, d)

uvTform = uvTformV;

if(size(d, 2) == 1)
    uvTform(:,7) = uvTformV(:,1)*d(1) + uvTformV(:,4)*d(2) + uvTformV(:,7);
    uvTform(:,8) = uvTformV(:,2)*d(1) + uvTformV(:,5)*d(2) + uvTformV(:,8);
    uvTform(:,9) = uvTformV(:,3)*d(1) + uvTformV(:,6)*d(2) + uvTformV(:,9);
else
    uvTform(:,7) = uvTformV(:,1).*d(:,1) + uvTformV(:,4).*d(:,2) + uvTformV(:,7);
    uvTform(:,8) = uvTformV(:,2).*d(:,1) + uvTformV(:,5).*d(:,2) + uvTformV(:,8);
    uvTform(:,9) = uvTformV(:,3).*d(:,1) + uvTformV(:,6).*d(:,2) + uvTformV(:,9);
end

% uvValidInd = uvTform(9,:) ~= 0;

uvTform = bsxfun(@rdivide, uvTform, uvTform(:, 9) + eps);
% uvTform(:,~uvValidInd) = 0;

end