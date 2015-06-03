function ifc=ifcvec(imorg,imdist);


% -----------COPYRIGHT NOTICE STARTS WITH THIS LINE------------
% Copyright (c) 2005 The University of Texas at Austin
% All rights reserved.
% 
% Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, 
% modify, and distribute this code (the source files) and its documentation for
% any purpose, provided that the copyright notice in its entirety appear in all copies of this code, and the 
% original source of this code, Laboratory for Image and Video Engineering (LIVE, http://live.ece.utexas.edu)
% at the University of Texas at Austin (UT Austin, 
% http://www.utexas.edu), is acknowledged in any publication that reports research using this code. The research
% is to be cited in the bibliography as:
% 
% H. R. Sheikh, A. C. Bovik, and G. de Veciana, "An Information Fidelity Criterion for Image 
% Quality Assessment Using Natural Scene Statistics," IEEE Transactios on Image Processing, in publication, May 2005.
% 
% 
% IN NO EVENT SHALL THE UNIVERSITY OF TEXAS AT AUSTIN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS DATABASE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF TEXAS
% AT AUSTIN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% THE UNIVERSITY OF TEXAS AT AUSTIN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE DATABASE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
% AND THE UNIVERSITY OF TEXAS AT AUSTIN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
% 
% -----------COPYRIGHT NOTICE ENDS WITH THIS LINE------------
%
%This is an implementation of the algorithm for calculating the
%Information Fidelity Criterion (IFC) between two images. Please refer
%to the following paper:
%
% H. R. Sheikh, A. C. Bovik, and G. de Veciana, "An Information Fidelity Criterion for Image 
% Quality Assessment Using Natural Scene Statistics," IEEE Transactios on Image Processing, in publication, May 2005.
%Download manuscript draft from http://live.ece.utexas.edu in the
%Publications link
%
%This implementation is slightly differnet from the one used to report
%results in the paper above. The modification have to do with using more
%subands than those used in the paper, better handling of image boundaries,
%and a window that automatically resizes itself based on the scale.
%
%Report bugfixes and comments to hamid.sheikh@ieee.org
%
%----------------------------------------------------------------------
% Prerequisites: The Steerable Pyramid toolbox. Available at
% http://www.cns.nyu.edu/~lcv/software.html
%
%Input : (1) img1: The reference image
%        (2) img2: The distorted image (order is important)
%
%Output: (1) The Information Fidelity between the two images

%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   ifc = ifcvec(img1, img2);
%
%Advanced Usage:
%   Users may want to modify the parameters in the code. 
%   (1) MxM is the block size that denotes the size of a vector used in the
%   GSM model.
%   (2) subbands included in the computation
%========================================================================

M=3;
subbands=[4 7 10 13 16 19 22 25];

% Do wavelet decomposition. This requires the Steerable Pyramid. You can
% use your own wavelet as long as the cell arrays org and dist contain
% corresponding subbands from the reference and the distorted images
% respectively.
[pyr,pind] = buildSpyr(imorg, 4, 'sp5Filters', 'reflect1'); % compute transform
org=ind2wtree(pyr,pind); % convert to cell array
[pyr,pind] = buildSpyr(imdist, 4, 'sp5Filters', 'reflect1');
dist=ind2wtree(pyr,pind);

% calculate the parameters of the distortion channel
[g_all,vv_all]=distsub_est_m(org,dist,subbands,M);

% calculate the parameters of the reference image
[ssarr, larr, cuarr]=refparams_vecgsm(org,subbands,M);

% reorder subbands. This is needed since the outputs of the above functions
% are not in the same order
vvtemp=cell(1,max(subbands));
ggtemp=vvtemp;
for(kk=1:length(subbands))
    vvtemp{subbands(kk)}=vv_all{kk};
    ggtemp{subbands(kk)}=g_all{kk};
end


% compute reference and distorted image information from each subband
for i=1:length(subbands)
    sub=subbands(i);
    g=ggtemp{sub};
    vv=vvtemp{sub};
    ss=ssarr{sub};
    lambda = larr(sub,:);, 
    cu=cuarr{sub};

    % how many eigenvalues to sum over. default is 1.
    neigvals=1;
    
    % compute the size of the window used in the distortion channel estimation, and use it to calculate the offset from subband borders
    % we do this to avoid all coefficients that may suffer from boundary
    % effects
    lev=ceil((sub-1)/6);
    winsize=2^lev+1; offset=(winsize-1)/2;
    offset=ceil(offset/M);
    
    
    % select only valid portion of the output.
    g=g(offset+1:end-offset,offset+1:end-offset);
    vv=vv(offset+1:end-offset,offset+1:end-offset);
    ss=ss(offset+1:end-offset,offset+1:end-offset);
    
    
    %IFC
    temp1=0; temp2=0;
    for j=1:length(lambda)
        temp1=temp1+sum(sum((log2(1+g.*g.*ss.*lambda(j)./(vv+1e-10))))); % IFC for the i'th subband. tolerence for zero variance
    end
    num(i)=temp1;
    
end

% compuate IFC and normalize to size of the image
ifc=sum(num)/prod(size(imorg));

