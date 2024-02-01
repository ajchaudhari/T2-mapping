<h1>Quantitative T2 mapping-based longitudinal assessment of brain injury and therapeutic rescue in the rat following acute organophosphate intoxication.<h1>

<h3>Code for T2 mapping using 2 methods:<h3>

<h4>1. Voxel-wise (curvefitting21.m)<h4>

<h4>2. Regional (curvefitting31.m)- averaging all voxels within a region and then computing T2<h4>

<h3>Note: MSME (multiple spin multiple echo) dataset is 4 dimensional- 9 slices across 15 timepoints (echo times/TEs). Segmentations were manually performed on T2-weighted scans<h3>

<h3>Below is the order of steps for the voxel-wise code (T2mapping_method2_mono.m):<h3>

<h4>1. Cropping_images.m is a function that performs cropping based on the center of mass algorithm<h4>

<h4>2. Labels are isolated and multiplied with MSME data using function Individual_ROI_mask.m<h4>

<h4>3. Monoexponential curve fitting  is then performed using function curvefitting21.m<h4>

<h4>4. Coefficients from the fit are then exported into an excel sheet<h4>

![image](https://github.com/aljesal/T2_mapping/assets/80182193/09a316d2-c676-4432-9f12-b7fcfa78615f)



