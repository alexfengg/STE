## MATLAB Code for CVPR '24 Paper: "A Subspace-Constrained Tyler’s Estimator and its Applications to Structure from Motion" 

This MATLAB code corresponds to the paper "A Subspace-Constrained Tyler’s Estimator and its Applications to Structure from Motion" by Feng Yu, Teng Zhang, and Gilad Lerman, presented at CVPR '24. The code contains two experiments: "Robust Fundamental Matrix Estimation" and "Initial Camera Removal for SfM".

### Photo Tourism Dataset

#### SfM Data:
The SfM data utilized in this study is sourced from the original Photo Tourism dataset provided by Soumyadip Sengupta [1]. The data is accessible through the "SfM_CVPR_code" repository. For detailed information, refer to the file "SfM_CVPR_code_data_readme.txt". The data can be downloaded from [this link](https://drive.google.com/file/d/1-nj44wLfFfZA8gu9fHxok6iLTGTbILUM/view?usp=sharing). After downloading, ensure that the data folder is located at "STE\CamRemoval_SfM\data\SfM_data".

#### SfM Clean Data:
To enhance program efficiency, the SfM data is processed and stored separately in another MATLAB file named "sfm_clean". Each entry in "sfm_clean" represents a pair of images (nodes). The data can be downloaded from [this link](https://drive.google.com/file/d/1-uOciycEpK04Sf_gxhs_qet1VivL05Xw/view?usp=sharing). After downloading, place the data folder at "STE\FundMatrix\data\sfm_clean". Each entry (datum{k}) includes the following information:
- Image pair indices: [i, j] = datum{k}.nodes
- Point correspondences: [x; x'] = datum{k}.points
- Calibration matrices: [K, K'] = datum{k}.K
- Relative rotation (R0): R0 = datum{k}.R0 (computed by Ri * Rj', where Ri and Rj are true rotations from SfM_data)
- Relative translation (t0): t0 = datum{k}.t0 (computed by Ri' * (tj - ti), where ti and tj are true translations from SfM_data)

### External Packages
1. Fast, Robust and Non-convex Subspace Recovery (FMS): [GitHub Repository](https://github.com/twmaunu/FMS)
2. A New Rank Constraint on Multi-view Fundamental Matrices, and its Application to Camera Location Recovery (SfM_CVPR2017_code): [Download Link](https://www.cs.unc.edu/~ronisen/SfM_CVPR2017_code.zip)
