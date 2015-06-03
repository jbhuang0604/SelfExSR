Oliver Woodford's Image-Based Rendering and Stereo algorithms, version 2.3.
----------------------------------------------------------------------

References
----------

This software package contains implementations of the image-based rendering algorithms described in the following papers:

[1] "Efficient New-view Synthesis using Pairwise Dictionary Priors", O. J. Woodford, I. D. Reid & A. W. Fitzgibbon. In proc. CVPR 2007.

[2] "On New View Synthesis Using Multiview Stereo", O. J. Woodford, I. D. Reid, P. H. S. Torr & A. W. Fitzgibbon. In proc. BMVC 2007.

[3] "Global Stereo Reconstruction under Second Order Smoothness Priors",  O. J. Woodford, P. H. S. Torr, I. D. Reid & A. W. Fitzgibbon. In proc. CVPR 2008.

License
-------

The software is licensed as described in the LICENSE.TXT file found in the root folder of this package.

Installing
----------

The software is Matlab source code and mex files. To install, unzip the files to a suitable location, maintaining the directory structure.

Initialise your Matlab path with the relevant directories for this software by calling startup.m, which resides in the root folder of this package.

Mex files will be compiled automatically when first called, but you must have a C++ compiler setup in Matlab for this to work. Do this by calling

   mex -setup

Preparing sequences
-------------------

This software only runs on calibrated image sequences - calibration software is not included. Examples of free and commercial calibration software include:
   Voodoo Camera Tracker
   Boujou

Each input sequence should be stored in its own directory, along with a data.mat file which contains both the projection matrices of the input views and either disparity values to sample at, or a set of 3-d scene points in world coordinates.

Projection matrices should be 3x4 matrices which convert homogenised world coordinates into image coordiantes of the form z*[x y 1]. They should encode the centre of the top left pixel of each input image as (0,0). The array of projection matrices, "Pi", should be 3x4xNumInputIms.

Disparity values, disps, should be a vector of positive, real-valued disparities at which to sample.

3-d points, points, should be a 3xNumPoints matrix of points, in world coordinates, which lie on objects within the scene. These points should be output by a structure from motion algorithm.

Matlab's current directory should then be set to that of the desired sequence prior to running the software.

Running
-------

Running this software is a four stage process:
1) Initialise the rendering parameters (including which method to use)
2) Initialise the projection matrices for the output views
3) Cache the input data for the current output view
4) Render/depth map the current output view

In detail:
1) To initialise parameters, call

   options = ojw_default_options(string);

where string is 'cvpr07' to use the algorithm from [1], 'bmvc07' to use that of [2], and 'cvpr08' to use that of [3].

The values in options can then be changed to tweak the algorithm.

2) To initialise output projection matrices, call

   Pout = ojw_genview(type, ind, offset);

where the input parameters are explained by calling

   help ojw_genview

Generally, all of the included functions return help text in this way.

3) Render/depth map the output view by calling

    [A out] = ibr_render(options, Pout);

These steps have been implemented in the example file ojw_script.m.

Bugs
----

If you find a bug in this software, excluding that software which has been authored outside the University of Oxford, please send a report of the bug to:

   ojw@robots.ox.ac.uk

Changes
-------

v2.0: Added implementation of stereo algorithm described in [3].
v2.1: Fixed a bug in vgg_interp2 which can cause crashes.
v2.2: Incorporate QPBO v1.1 source (previously v1.0).
v2.3: Added function to automatically download Middlebury sequences. Changed setup process somewhat.
v2.4: Added truncated quadratic kernel to stereo code. Fixed bugs kindly reported by Dhruv Batra and Petter Strandmark.

