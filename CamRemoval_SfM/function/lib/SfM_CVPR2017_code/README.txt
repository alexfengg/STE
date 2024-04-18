%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%

SfM_CVPR_code/Data : contains .mat files for each of the 14 scenes.
Note : there are many duplicates and not necessary variables

Important variables :

n - number of images
t,R - Ground Truth(GT) camera location and orientation
keep/AdjMat - logical matrix of size nxn, keep(i,j)=1 means pair-(i,j) is estimated.
G_gt - GT pairwise Rotation (R_i x R_j')
E/E_gt - collection of GT essential matrices = 
H_mat - Estimated pair-wise Rotation R_ij=(R_i x R_j')
tijGT - Estimated pair-wise translation t_ij
E_est -Estimated pair-wise Essentials [t_ij]_{x} R_ij
corr - contains corresponding points for every pair of image. corr{i,j}(:,:,1) for img-i and corr{i,j}(:,:,2) for img-j.
K - Camera Calibration

For using any new data you need to make sure the format matches with the code.

%%%%%%%%%%%% Dependency %%%%%%%%%%%%%%
1) Install CVX and make sure the path is accessible. You can add : ‘addpath(genpath(PATH_TO_CVX))’ at the beginning of script_run_r3.m
2) I provide only 1 processed data file ‘Alamo’ along with the code. Rest of the processed data files can be downloaded separately. You need to place each of those .mat files in ‘Data’ folder.


%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%
Main file : script_run_r3.m (See explanation in code)

Try : script_run_r3(1,1,50);

%%%%%%%% License %%%%%%%%%%%%%%%
Feel free to use the code for academic research purposes only.
Please DONOT re-distribute.

%%%%%%%% Acknowledgement %%%%%%%%%%
The code involves a pipeline of Structure from Motion built by several authors.
(A) The data is originally collected by the authors in [1]. We just read, process and store the data as .mat file.
(B) The code for the Baseline algorithm (LUD) is developed by Onur Ozyesil [2].
(C) The code for IRLA-ADMM based rank-contsrained optimization on Fundamental matrices is developed by Soumyadip Sengupta [3].

%%%%%%%% Contact %%%%%%%%
For any query, contact : 
Soumyadip Sengupta
University of Maryland, College Park
sengupta@umd.edu.

%%%%%%% Reference %%%%%%
[1] Kyle Wilson and Noah Snavely. "Robust Global Translations with 1DSfM." ECCV 2014.

[2] Onur Ozyesil and Amit Singer. "Robust camera location estimation by convex programming." CVPR 2015.

[3] S Sengupta, T Amir, M Galun, T Goldstein, DW Jacobs, A Singer, R Basri. "A New Rank Constraint on Multi-view Fundamental Matrices, and its Application to Camera Location Recovery." CVPR 2017
