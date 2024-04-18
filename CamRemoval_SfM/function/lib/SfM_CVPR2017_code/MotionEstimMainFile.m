%% This is the main file demonstrating the use of our algorithms for the
%% paper titled "Robust Camera Location Estimation by Convex Programming".
%%
%% IMPORTANT: These codes are not pre-tested, hence, DO NOT DISTRIBUTE...
%%
%% These codes are a part of private communication between Onur Ozyesil 
%% and Soumyadip Sengupta.

%% Load your data: Please note that a specific data structure is assumed.
%%
%% However, be aware of the fact that the CHOICE OF GLOBAL COORDINATE 
%% SYSTEMS, MEASUREMENT COORDINATE SYSTEM (e.g., center of each image, x&y 
%% directions, etc.), TYPE OF MEASUREMENTS (e.g., one may be assuming 
%% R_ij = R_i^(-1)*R_j or R_ij = R_i*R_j^(-1) or R_ij = R_j*R_i^(-1), etc.)
%% may be (and in many cases are) different for each dataset, which may 
%% produce meaningless estimates, comparison results with our methods. 
%% Please refer to our paper "Robust Camera Location Estimation by Convex 
%% Programming" for our measurement types, corrdinate system choices, etc.,
%% and also refer to the files for the datasets you will be using for 
%% comparison.

%% FILES ARE WRITTEN AND (INSUFFICIENTLY) TESTED USING MATLAB R2012a. 
%% Some function may need to be updated for newer versions of MATLAB.

clear all; 
clc;

Data = load('PATHOFYOURDATAFILE');

%% Prepare data
SO3Mats_orig = Data.CoreData.SO3Mats_orig; % Ground truth camera orientations Ri
t_orig = Data.CoreData.t_orig; % Ground truth camera locations
K_Mats = Data.CoreData.K_Mats; % Camera calibration matrices
CorrPts = Data.CoreData.CorrPnts; % Corresponding feature point extractions
AdjMat = Data.CoreData.AdjMat; % Adjacency matrix of corresponding images, we assume AdjMat is CONNECTED
Hmat = Data.CoreData.H; % Large matrix of pairwise rotations, ij'th 3by3 block of H is the estimate of Ri'*Rj
n = size(AdjMat,1);

%% Orientation estimation by iterative SO3 averaging:
tRotEst_s = tic;
iterNum_IEVM = 5; % This choice, i.e. iterNum_IEVM = 5, is completely 
                  % arbitrary, make your own "heuristic" choice
[SO3Mats_est, IndVec_IAV] = IterativeSO3Average(AdjMat,Hmat,iterNum_IEVM);
% SO3Mats_est = permute(SO3Mats_est,[2 1 3]);
tRotEst = toc(tRotEst_s);

%% Evaluate orientation estimation performance for camera indices returned by IEVM
[SO3Mats_corrected, MSE_Error_EVM, R_global] = GlobalSOdCorrectLeft(SO3Mats_est, SO3Mats_orig(:,:,IndVec_IAV));


%% Cut the data using IndVec_IEVM
IEVM_cut = zeros(n,1); IEVM_cut(IndVec_IAV) = 1; IEVM_cut = logical(IEVM_cut);
n = sum(IEVM_cut);
AdjMat = AdjMat(IEVM_cut,IEVM_cut);
Hmat = Hmat(logical(kron(IEVM_cut,ones(3,1))),logical(kron(IEVM_cut,ones(3,1))));
CorrPts = CorrPts(IEVM_cut,IEVM_cut);
K_Mats = K_Mats(:,:,IEVM_cut);
SO3Mats_orig = SO3Mats_orig(:,:,IEVM_cut);
t_orig = t_orig(:,IEVM_cut);

%% Check for rigidity:
tParRig_s = tic;
rigidityDec = ParallelRigidityTest(AdjMat,3);
if rigidityDec
    fprintf('\n Measurement graph is parallel rigid :) \n');
    PRcomp = logical(ones(n,1));
else
    fprintf('\n Measurement graph is not parallel rigid! \n Extracting largest maximally parallel rigid component... \n');
    [LarMaxRigCompInds,~] = LargestMaxParRigComp(AdjMat,3);
    PRcomp = logical(zeros(n,1)); PRcomp(LarMaxRigCompInds) = true;
    n = sum(PRcomp);
    AdjMat = AdjMat(PRcomp,PRcomp);
    Hmat = Hmat(logical(kron(PRcomp,ones(3,1))),logical(kron(PRcomp,ones(3,1))));
    CorrPts = CorrPts(PRcomp,PRcomp);
    SO3Mats_orig = SO3Mats_orig(:,:,PRcomp);
    t_orig = t_orig(:,PRcomp);
    SO3Mats_est = SO3Mats_est(:,:,PRcomp);
end
tParRig = toc(tParRig_s);

%% Estimate pairwise directions:
invK_Mats = zeros(3,3,n);
for k = 1:n
    invK_Mats(:,:,k) = [1/K_Mats(1,1,k) 0 -K_Mats(1,3,k)/K_Mats(1,1,k);...
                    0 1/K_Mats(2,2,k) -K_Mats(2,3,k)/K_Mats(2,2,k); 0 0 1];
end

tRobDir_s = tic;
tijMat = RobustSignedLineEstimInvK(SO3Mats_est,CorrPts,invK_Mats);
tRobDir = toc(tRobDir_s);

%% Estimate locations by LUD
% Use IRLS:
optsIRLS.tolIRLS = 1e-5;       
optsIRLS.tolQuad = 1e-10;
optsIRLS.delt = 1e-16;
optsIRLS.maxit   = 200;  
[t_LUD, out_LUD] = LocationEstimByLUD(AdjMat,tijMat,optsIRLS);

% Or use SOCP solver of the CVX package: 
% !!! NOTE THAT THIS CHOICE REQUIRES THE CVX PACKAGE TO
% BE INSTALLED !!!
% [t_LUD, ~] = SynchDirectionsByLUDCvX(AdjMat, tijMat);
% %% If desired, locally refine
% p_LUD = 1.5;
% [t_LUD,~] = RefineLUDByFminuncon(vec(t_LUD),AdjMat,tijMat,p_LUD);

%% Evaluate performance:
[t_fit_LUD,~,~,NRMSE_LUD,~] = SimpleTransScaleRemove(R_global*reshape(t_LUD,3,n), t_orig,'L2');
figure; plot3(t_fit_LUD(1,:),t_fit_LUD(2,:),t_fit_LUD(3,:),'b*');
hold on; plot3(t_orig(1,:),t_orig(2,:),t_orig(3,:),'r*');


