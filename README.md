# Bone & bone lesion segmentation in SPECT/CT with the modified Robust Fuzzy C-means
This is an implementation of my paper:

<a href="https://www.spiedigitallibrary.org/conference-proceedings-of-spie/10949/109491W/Incorporating-CT-prior-information-in-the-robust-fuzzy-C-means/10.1117/12.2506805.short">Chen, Junyu, et al. "Incorporating CT prior information in the robust fuzzy C-means algorithm for QSPECT image segmentation." Proceedings Volume 10949, SPIE Medical Imaging 2019: Image Processing.</a>

This segmentation method was based the <a href="https://www.sciencedirect.com/science/article/pii/S1077314201909518">Robust Fuzzy C-means algorithm</a> with the incorporation of prior knowledge from the CT clustering.

## Preprocess Images:
Crop a small region of interest in both SPECT and CT image.

![](https://github.com/junyuchen245/SPECT-CT-Seg-RFCM/blob/master/sample_img/cropping_imgs.png)
## Objective Function:
![](https://github.com/junyuchen245/SPECT-CT-Seg-RFCM/blob/master/sample_img/objective_func.png)
## Sample Results:
![](https://github.com/junyuchen245/SPECT-CT-Seg-RFCM/blob/master/sample_img/seg_results.png)

Red contours in the images represent segmented bone regions and green contours indicate
lesion regions. The segmentation result in image (a) was obtained with β = 0.0005 and γ = 0. The segmentation results in images (b), (c), and (d) were obtained with β = 0.0005 and γ = 0.005, 0.01, and 0.05, respectively.
