function [ssarr, l_arr, cu_arr]=refparams_vecgsm(org,subands,M)

%This function computes the parameters of the reference image. This is
%called by vifvec.m.

for i=1:length(subands);
    sub=subands(i);
    y=org{sub};
    
    sizey=floor(size(y)./M)*M; % crop  to exact multiple size
    y=y(1:sizey(1),1:sizey(2));
    
    
    % Collect MxM blocks. Rearrange each block into an
    % M^2 dimensional vector and collect all such vectors.
    % Collece ALL possible MXM blocks (even those overlapping) from the subband
    temp=[];
    for j=1:M
        for k=1:M
            temp=cat(1,temp,reshape(y(k:end-(M-k), j:end-(M-j)),1,[]));
        end
    end
    
    % estimate mean and covariance
    mcu=mean(temp')';
    cu=((temp-repmat(mcu,1,size(temp,2)))*(temp-repmat(mcu,1,size(temp,2)))')./size(temp,2); % covariance matrix for U
    
    % Collect MxM blocks as above. Use ONLY non-overlapping blocks to
    % calculate the S field
    temp=[];
    for j=1:M
        for k=1:M
            temp=cat(1,temp,reshape(y(k:M:end, j:M:end),1,[]));
        end
    end

    % Calculate the S field
    ss=(inv(cu)*temp);
    ss=sum(ss.*temp)./(M*M);
    ss=reshape(ss,sizey/M);
    
    % Eigen-decomposition
    [v,d]=eig(cu);
    l_arr(sub,:)=diag(d)';
    
    % rearrange for output
    ssarr{sub}=ss;
    temp=0;
    d=diag(d);
    cu_arr{sub}=cu;
end    

