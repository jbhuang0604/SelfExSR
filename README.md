## Single Image Super-Resolution from Transformed Self-Exemplars (CVPR 2015)

### Introduction

This is the research code for the paper:

[Jia-Bin Huang](https://sites.google.com/site/jbhuang0604/), [Abhishek Singh](https://sites.google.com/site/abhishek486/), and [Narendra Ahuja] (http://vision.ai.illinois.edu/ahuja.html), "Single Image Super-Resolution from Transformed Self-Exemplars", CVPR 2015 [PDF](https://uofi.box.com/shared/static/8llt4ijgc39n3t7ftllx7fpaaqi3yau0.pdf)

The proposed algorithm achieves the state-of-the-art performance on image super-resolution *without* the need of any external training dataset, feature extraction and complicated learning algorithms. For more details, please visit our [Project page](https://sites.google.com/site/jbhuang0604/publications/struct_sr).

All the datasets (Set5, Set14, Urban 100, BSD 100, Sun-Hays 80), precomputed results and visual comparisons can be found in the following sections.

### Citation

If you find the code and dataset useful in your research, please consider citing:

    @inproceedings{Huang-CVPR-2015,
        title={Single Image Super-Resolution From Transformed Self-Exemplars},
        Author = {Huang, Jia-Bin and Singh, Abhishek and Ahuja, Narendra},
        booktitle = {Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition},
        pages={5197--5206},
        Year = {2015}
    }

### Contents
|  Folder    | description |
| ---|---|
|cache | cached data for vanishing point detection|
|data|Testing images of five datasets (Set5, Set14, Urban 100, BSD 100, Sun-Hays 80). All the images have been cropped according to the desired super-resolution factor. This avoids misalignment of the groundtruth high-resolution images and the super-resolved images|
|external|We use the vgg_interp2 from `imrender` to perform bilinear interpolation|
|quant_eval|Quantitative evaluation code|
|reference| A copy of the CVPR paper and the bibtex|
|source|MATLAB source code|

To run the algorithm on all datasets, simply run the `sr_demo_bacth.m`. Note that it is an educational code that is not optimized for speed. If timing is a concern, you can achieve visually similar results with small numbers of iterations, e.g., set the number of iterations `opt.numIter = 5;` in the file `sr_init_opt.m`. An example of the speed and quality trade-off can be found in Fig. 10 in the paper.

Feedbacks and comments are welcome! Feel free to contact me via jbhuang1@illinois.edu.

Enjoy!

Note: For vanishing point detection only windows executable is provided (from [Image Completion using Planar Structure Guidance](https://github.com/jbhuang0604/StructCompletion)), a cross-platform version will be included later.

### Comparison with the state-of-the-art

##### Datasets

The full super-resolution results on Set 5, Set 14, Urban 100, BSD 100 and Sun-Hays 80 are available.

| Dataset | Image source | Download full results |
|---- | ---|----| ----|
| **Set 5** |  [Bevilacqua et al. BMVC 2012](http://people.rennes.inria.fr/Aline.Roumy/results/SR_BMVC12.html)  | [link](https://uofi.box.com/shared/static/kfahv87nfe8ax910l85dksyl2q212voc.zip) (16.1 MB)|
| **Set 14** |  [Zeyde et al. LNCS 2010](https://sites.google.com/site/romanzeyde/research-interests)  | [link](https://uofi.box.com/shared/static/igsnfieh4lz68l926l8xbklwsnnk8we9.zip) (86.0 MB)|
| **Urban 100** | [Huang et al. CVPR 2015](https://sites.google.com/site/jbhuang0604/publications/struct_sr)  | [link](https://uofi.box.com/shared/static/65upg43jjd0a4cwsiqgl6o6ixube6klm.zip) (1.14 GB)|
| **BSD 100** | [Martin et al. ICCV 2001](https://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/) | [link](https://uofi.box.com/shared/static/qgctsplb8txrksm9to9x01zfa4m61ngq.zip) (568 MB)|
| **Sun-Hays 80** | [Sun and Hays ICCP 2012](http://cs.brown.edu/~lbsun/SRproj2012/SR_iccp2012.html) | [link](https://uofi.box.com/shared/static/rirohj4773jl7ef752r330rtqw23djt8.zip) (311 MB)|

Set 5 dataset - [link](https://uofi.box.com/shared/static/kfahv87nfe8ax910l85dksyl2q212voc.zip)
![Set 5](https://uofi.box.com/shared/static/sk9duzzu63x80zdgwszf1vwqwl0ea7zx.jpg)

Set 14 dataset - [link](https://uofi.box.com/shared/static/igsnfieh4lz68l926l8xbklwsnnk8we9.zip)
![Set 14](https://uofi.box.com/shared/static/b8067imlbojcdk6guepudlj7c0wh3kmd.jpg)

Urban 100 dataset - [link](https://uofi.box.com/shared/static/65upg43jjd0a4cwsiqgl6o6ixube6klm.zip)
![Urban 100](https://uofi.box.com/shared/static/20cy9kji3990py2jwu4uwidho3wh2ke0.jpg)

BSD 100 dataset - [link](https://uofi.box.com/shared/static/qgctsplb8txrksm9to9x01zfa4m61ngq.zip)
![BSD 100](https://uofi.box.com/shared/static/yx1eqfb2yewy5fj2bvnxxb766irgpfh4.jpg)

Sun-Hays 80 dataset - [link](https://uofi.box.com/shared/static/rirohj4773jl7ef752r330rtqw23djt8.zip)
![Sun-Hays 80](https://uofi.box.com/shared/static/5mal435jvm5tanszrd95e1orltomq8s3.jpg)

##### State-of-the-art image super-resolution algorithms

In each dataset, we include results of the state-of-the-art single image super-resolution algorithms:

| Image | Description|
| ---|---|
| **HR**| High-resolution images. All images were cropped so that each dimension is a multiplication of the super-resolution factor. This avoids the misalignment problem in the quantitative comparison. |
| **LR** | Low-resolution test images generated with bicubic kernel downsampling.| 
|**bicubic**|Bicubic interpolation|
|**nearest**|Nearest-neighbor interpolation|
|**SelfExSR**|Our result|
|**A+**|R. Timofte, V. De Smet, and L. Van Gool, A+: Adjusted Anchored Neighborhood Regression for Fast Super-Resolution. In Asian Conference on Computer Vision (ACCV 2014). Code available [here](http://www.vision.ee.ethz.ch/~timofter/ACCV2014_ID820_SUPPLEMENTARY/index.html)|
|**Abhishek**|A. Singh and N. Ahuja, Super-Resolution Using Sub-Band Self-Similarity. In Asian Conference on Computer Vision (ACCV 2014). No publicly implementation available. Results were provided by the authors.|
|**Kim**|K. I. Kim and Y. Kwon, “Single-image super-resolution using sparse regression and natural image prior”, IEEE Trans. Pattern Analysis and Machine Intelligence, vol. 32, no. 6, pp. 1127-1133, 2010. Code available [here](https://people.mpi-inf.mpg.de/~kkim/supres/supres.htm)|
|**Glasner**| Daniel Glasner, Shai Bagon, Michal Irani. Super-Resolution From a Single Image, In International Conference on Computer Vision (ICCV 2009). No public implementation available. Results were generated by our own implementation.|
|**ScSR**|Jianchao Yang, John Wright, Thomas Huang, and Yi Ma. Image super-resolution via sparse representation. IEEE Transactions on Image Processing, Vol 19, Issue 11, pp2861-2873, 2010. Code available [here](http://www.ifp.illinois.edu/~jyang29/ScSR.htm)|
|**SRCNN**| Chao Dong, Chen Change Loy, Kaiming He, Xiaoou Tang. Learning a Deep Convolutional Network for Image Super-Resolution, in European Conference on Computer Vision (ECCV 2014). Code available [here](http://mmlab.ie.cuhk.edu.hk/projects/SRCNN.html)|

### Qualitative comparison

In our supplementary material, we includde 120 sample comparisons with the state-of-the-art algorithms. Download the document [here](https://uofi.box.com/shared/static/k0a2wziaavwoo557aygkbwsa6syln1b9.pdf).

You can browse and compare our results with other methods via the following links.
 * [Urbana 100](https://dl.dropboxusercontent.com/u/2810224/Homepage/publications/2015/SuperResolution_CVPR_2015/supp/Urban_SRF_4.html) - Super-resolution factor 4x
 * [BSD 100](https://dl.dropboxusercontent.com/u/2810224/Homepage/publications/2015/SuperResolution_CVPR_2015/supp/BSD_SRF_3.html) - Super-resolution factor 3x
 * [Sun-Hays 80](https://dl.dropboxusercontent.com/u/2810224/Homepage/publications/2015/SuperResolution_CVPR_2015/supp/Sun_Hays_SRF_8.html) - Super-resolution factor 8x

### Quantitative comparisons

We report three types of metrics
  * PSNR: [Peak signal-to-noise ratio](http://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio)
  * SSIM: [Structural similarity index](https://ece.uwaterloo.ca/~z70wang/research/ssim/)
  * IFC:  [Information fidelity criterion](http://live.ece.utexas.edu/research/quality/)

##### Results on Set 5

|  Scale    | Bicubic | ScSR  | Kim | Sub-band |  Glasner |SRCNN  | A+ | Ours |
|:---------:|:-------:|:--------:|:------:|:------------:|:---------:|:--------:|:------:|:----:|
| **2x** - PSNR|   33.64	|   35.78	|   36.24	|    Sub-band	|   35.43	|   36.28	|    A+	|   36.50	| 
| **3x** - PSNR|   30.39	|   31.34	|   32.30	|    Sub-band	|   31.10	|   32.37	|    A+	|   32.62	|
| **4x** - PSNR|   28.42	|   29.07	|   30.07	|    Sub-band	|   28.84	|   30.08	|    A+	|   30.33	|
||
| **2x** - SSIM|  0.9292	|  0.9485	|  0.9518	|  Sub-band	|  0.9452	|  0.9509	|  A+	|  0.9537	| 
| **3x** - SSIM|  0.8678	|  0.8869	|  0.9041	|  Sub-band	|  0.8811	|  0.9025	|  A+	|  0.9094	| 
| **4x** - SSIM|  0.8101	|  0.8263	|  0.8553	|  Sub-band	|  0.8210	|  0.8525	|  A+	|  0.8623	| 
||
| **2x** - IFC|    5.72	|    6.94	|    7.05	|    Sub-band	|    6.70	|    6.85	|    A+	|    7.83	| 
| **3x** - IFC|    3.45	|    3.98	|    4.25	|    Sub-band	|    3.68	|    4.11	|    A+	|    4.76	| 
| **4x** - IFC|    2.28	|    2.57	|    2.82	|    Sub-band	|    2.42	|    2.76	|    A+	|    3.19	| 

##### Results on Set 14 

|  Scale    | Bicubic | ScSR  | Kim | Sub-band |  Glasner |SRCNN  | A+ | Ours |
|:---------:|:-------:|:--------:|:------:|:------------:|:---------:|:--------:|:------:|:----:|
| **2x** - PSNR|   30.22	|   31.64	|   32.14	|    Sub-band	|   31.41	|   32.00	|    A+	|   32.23	| 
| **3x** - PSNR|   27.53	|   28.19	|   28.96	|    Sub-band	|   28.21	|   28.90	|    A+	|   29.16	| 
| **4x** - PSNR|   25.99	|   26.40	|   27.18	|    Sub-band	|   26.43	|   27.13	|    A+	|   27.40	| 
||
| **2x** - SSIM|  0.8683	|  0.8940	|  0.9031	|  Sub-band	|  0.8881	|  0.9012	|  A+	|  0.9036	| 
| **3x** - SSIM|  0.7737	|  0.7977	|  0.8140	|  Sub-band	|  0.7926	|  0.8124	|  A+	|  0.8197	| 
| **4x** - SSIM|  0.7023	|  0.7218	|  0.7434	|  Sub-band	|  0.7163	|  0.7395	|  A+	|  0.7518	| 
||
| **2x** - IFC|    5.74	|    6.83	|    6.92	|    Sub-band	|    6.47	|    6.68	|    A+	|    7.60	| 
| **3x** - IFC|    3.33	|    3.75	|    3.92	|    Sub-band	|    3.59	|    3.81	|    A+	|    4.38	| 
| **4x** - IFC|    2.18	|    2.46	|    2.57	|    Sub-band	|    2.30	|    2.50	|    A+	|    2.90	| 

##### Results on Urban 100

|  Scale    | Bicubic | ScSR  | Kim | Sub-band |  Glasner |SRCNN  | A+ | Ours |
|:---------:|:-------:|:--------:|:------:|:------------:|:---------:|:--------:|:------:|:----:|
| **2x** - PSNR|   26.66	|   28.26	|   28.74	|   28.34	|   27.85	|   28.65	|   28.87	|   29.38	| 
| **4x** - PSNR|   23.14	|   24.02	|   24.20	|   24.19	|   23.58	|   24.14	|   24.34	|   24.82	| 
||
| **2x** - SSIM|  0.8408	|  0.8828	|  0.8940	|  0.8820	|  0.8709	|  0.8909	|  0.8957	|  0.9032	| 
| **4x** - SSIM|  0.6573	|  0.7024	|  0.7104	|  0.7115	|  0.6736	|  0.7047	|  0.7195	|  0.7386	|
||
| **2x** - IFC|    5.72	|    6.98	|    6.86	|    7.08	|    6.17	|    6.66	|    8.02	|    7.96	| 
| **4x** - IFC|    2.27	|    2.75	|    2.71	|    2.72	|    2.35	|    2.63	|    3.16	|    3.33	|

##### Results on BSD 100

|  Scale    | Bicubic | ScSR  | Kim | Sub-band |  Glasner |SRCNN  | A+ | Ours |
|:---------:|:-------:|:--------:|:------:|:------------:|:---------:|:--------:|:------:|:----:|
| **2x** - PSNR|   29.55	|   30.77	|   31.11	|   30.73	|   30.28	|   31.11	|   31.22	|   31.18	| 
| **3x** - PSNR|   27.20	|   27.72	|   28.17	|   27.88	|   27.06	|   28.20	|   28.30	|   28.30	| 
| **4x** - PSNR|   25.96	|   26.61	|   26.71	|   26.60	|   26.17	|   26.70	|   26.82	|   26.85	| 
||
| **2x** - SSIM|  0.8425	|  0.8744	|  0.8840	|  0.8774	|  0.8621	|  0.8835	|  0.8862	|  0.8855	| 
| **3x** - SSIM|  0.7382	|  0.7647	|  0.7788	|  0.7714	|  0.7368	|  0.7794	|  0.7836	|  0.7843	| 
| **4x** - SSIM|  0.6672	|  0.6983	|  0.7027	|  0.7021	|  0.6747	|  0.7018	|  0.7089	|  0.7108	| 
||
| **2x** - IFC|    5.26	|    6.20	|    6.30	|    6.36	|    5.56	|    6.09	|    7.15	|    6.84	| 
| **3x** - IFC|    3.00	|    3.37	|    3.49	|    3.17	|    2.72	|    3.39	|    3.92	|    3.81	| 
| **4x** - IFC|    1.91	|    2.22	|    2.20	|    2.18	|    1.86	|    2.18	|    2.51	|    2.46	| 
