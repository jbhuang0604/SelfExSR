function [g_all, vv_all]=vifsub_est_M(org,dist, subbands, M); 

% uses convolution for determining the parameters of the distortion channel
% Called by vifvec.m

tol = 1e-15; % tolernace for zero variance. Variance below this is set to zero, and zero is set to this value to avoid numerical issues.


for i=1:length(subbands)
    sub=subbands(i);
    y=org{sub};
    yn=dist{sub};

    % compute the size of the window used in the distortion channel estimation
    lev=ceil((sub-1)/6);
    winsize=2^lev+1; offset=(winsize-1)/2;
    win = ones(winsize);
    
    % force subband size to be multiple of M
    newsize=floor(size(y)./M)*M;
    y=y(1:newsize(1),1:newsize(2));
    yn=yn(1:newsize(1),1:newsize(2));

    % Correlation with downsampling. This is faster than downsampling after
    % computing full correlation.
    winstep=[M M];
    winstart=[1 1].*floor(M/2)+1;
    winstop=size(y)-ceil(M/2)+1;
    
    % mean
    mean_x = corrDn(y,win/sum(win(:)),'reflect1',winstep, winstart,winstop);
    mean_y = corrDn(yn,win/sum(win(:)),'reflect1',winstep, winstart,winstop);
    % cov
    cov_xy = corrDn(y.*yn, win, 'reflect1',winstep, winstart,winstop) - sum(win(:)).*mean_x.*mean_y;
    % var
    ss_x = corrDn(y.^2,win, 'reflect1',winstep, winstart,winstop) - sum(win(:)).*mean_x.^2;
    ss_y = corrDn(yn.^2,win, 'reflect1',winstep, winstart,winstop) - sum(win(:)).*mean_y.^2;

    
    % get rid of numerical problems, very small negative numbers, or very
    % small positive numbers, or other theoretical impossibilities.
    ss_x(ss_x<0)=0;
    ss_y(ss_y<0)=0;
   
    % Regression 
    g = cov_xy./(ss_x+tol);
    
    % Variance of error in regression
    vv = (ss_y - g.*cov_xy)/(sum(win(:)));
    
    % get rid of numerical problems, very small negative numbers, or very
    % small positive numbers, or other theoretical impossibilities.
    g (ss_x < tol) = 0;
    vv (ss_x < tol) = ss_y (ss_x < tol);
    ss_x(ss_x<tol)=0;
    
    g (ss_y < tol) = 0;
    vv (ss_y < tol) = 0;
    
    % constrain g to be non-negative. 
    vv(g<0)=ss_y(g<0);
    g(g<0)=0;
    
    % take care of numerical errors, vv could be very small negative
    vv( vv <= tol) = tol;
    
    g_all{i}=g;
    vv_all{i}=vv;
    
end