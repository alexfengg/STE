## MATLAB Code for CVPR '24 Paper: "A Subspace-Constrained Tyler’s Estimator and its Applications to Structure from Motion" 

MATLAB code for the paper "A Subspace-Constrained Tyler’s Estimator and its Applications to Structure from Motion" [1]. Two experiments, "Robust Fundamental Matrix Estimation (RFME)" and "Initial Camera Removal for SfM (ICRS)", are provided, corresponding Section 4.1 and Section 4.2 in paper. 

The experiments are based on the [Photo Tourism Dataset](https://www.cs.cornell.edu/projects/1dsfm/). The SfM data (utilized in the ICRS) was generated by [2] and one of the data set is accessible through the ["SfM_CVPR_code"](https://www.cs.unc.edu/~ronisen/SfM_CVPR2017_code.zip) repository. For detailed information, refer to the file "SfM_CVPR_code\README.txt". 

#### Experiment 1: Robust Fundamental Matrix Estimation

##### Step 1: download the SfM (Clean) Data
Download the processed SfM data from [this link](https://drive.google.com/file/d/1-nj44wLfFfZA8gu9fHxok6iLTGTbILUM/view?usp=sharing) and place the "sfm_clean" folder at "STE\FundMatrix\data\sfm_clean".

##### Step 2: implementation
Run the MATLAB script run_all.m located in the directory "STE\FundMatrix" to produce the errors of the Subspace-Constrained Tyler's Estimator (STE) and other Robust Subspace Recovery (RSR) methods of TME, FMS, and SFMS.

```
run_all
```
Generate the figures by calling the script 'plotting.m'.
```
plotting
```

#### Experiment 2: Initial Camera Removal for SfM

##### Step 1: download the  SfM Data
Download the required SfM data from [this link](https://drive.google.com/file/d/1-uOciycEpK04Sf_gxhs_qet1VivL05Xw/view?usp=sharing) and place the "SfM_data" folder at "STE\CamRemoval_SfM\data\SfM_data".

##### Step 2: implementation
Follow the instructions below to produce the errors of STE and other RSR methods (TME, FMS, SFMS):
1. Run the LUD pipeline [6] by calling the function run_baseline located in the directory "STE\CamRemoval_SfM".
```
run_baseline
```
2. Generate the stacked estimated essential matrix and fill the zero blocks by running the matrix completion algorithm.
```
run_MatComp
```
3. Run the LUD+RSR pipeline by calling the function run_all located in the directory "STE\CamRemoval_SfM".
```
run_all
```

### Acknowledgement
1. The data is originally collected by the authors in [3]. The data are processed and stored as .mat file. 
2. SfM data was processed by the authors in [2] and it was further and stored separately in "sfm_clean". Note that only correspondence, calibration and the (ground-truth) relative pose are needed in estimating the fundamental matrix. Only storing these data will greatly increase the computational efficiency. Each entry in "sfm_clean" represents a pair of images (nodes). 
Note: Each entry (datum{k}) includes the following information:
- Image pair indices: [i, j] = datum{k}.nodes
- Point correspondences: [x; x'] = datum{k}.points
- Calibration matrices: [K, K'] = datum{k}.K
- Relative rotation (R0): R0 = datum{k}.R0 (computed by Ri * Rj', where Ri and Rj are true rotations from SfM_data)
- Relative translation (t0): t0 = datum{k}.t0 (computed by Ri' * (tj - ti), where ti and tj are true translations from SfM_data)
3. Parallel Computing: it is recommended that run the parfor in thread-based enviroments for experiment 1. To setup, change the setting in 'Preferences/Parallel Computing Toolbox/Parallel Environment/Default Profile' from "Processes" to "Threads".
The difference between the process-based (default in Matlab) environments and thread-based environments see [here](https://www.mathworks.com/help/parallel-computing/choose-between-thread-based-and-process-based-environments.html#mw_6bbf0761-74c0-404e-9db6-77b82c7c138c). 
4. The comparing robust subspace recovery (RSR) methods include FMS [4] and TME [5].



### References
[1] Feng Yu, Teng Zhang, and Gilad Lerman. "A Subspace-Constrained Tyler's Estimator and its Applications to Structure from Motion." arXiv:2404.11590.

[2] Soumyadip Sengupta, Tal Amir, Meirav Galun, Tom Goldstein, David Jacobs, Amit Singer, and Ronen Basri. "A new rank constraint on multi-view fundamental matrices, and its application to camera location recovery." CVPR 2017.

[3] Kyle Wilson and Noah Snavely. "Robust Global Translations with 1DSfM." ECCV 2014. 

[4]Gilad Lerman and Tyler Maunu. "Fast, Robust and Non-convex Subspace Recovery." Information and Inference: A Journal of the IMA 2018.

[5] Teng Zhang. "Robust subspace recovery by Tyler’s M-estimator." Information and Inference: A Journal of the IMA 2016.

[6] Onur Ozyesil and Amit Singer. "Robust camera location estimation by convex programming." CVPR 2015.