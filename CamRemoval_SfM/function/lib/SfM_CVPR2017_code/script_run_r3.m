function []=script_run_r3(id,rand_id,N)
%% Runs Baseline (LUD) and Our method and compares both of them
%Input : 
%id : id of dataset
%N : Number of images to use
%rand_id : random number generator for selecting different sub-samples of N
%images. To Duplicate the result of the paper. use rand_id 1:5 inside
%matlab parpool environment. The rng(rand_id) generates different sequence
%in normal and parallel processing environment.

MAXIT=20; %Number of iteration of IRLS
%% Load data
DataSetName1{1} = 'Alamo';
DataSetName1{2} = 'Ellis_Island';
DataSetName1{3} = 'Madrid_Metropolis';
DataSetName1{4} = 'Montreal_Notre_Dame';
DataSetName1{5} = 'Notre_Dame';
DataSetName1{6} = 'NYC_Library';
DataSetName1{7} = 'Piazza_del_Popolo';
DataSetName1{8} = 'Piccadilly';
DataSetName1{9} = 'Roman_Forum';
DataSetName1{10} = 'Tower_of_London';
DataSetName1{11} = 'Trafalgar';
DataSetName1{12} = 'Union_Square';
DataSetName1{13} = 'Vienna_Cathedral';
DataSetName1{14} = 'Yorkminster';
DataSetName1{15} = 'Gendarmenmarkt';

DataSetName=DataSetName1{id};

load(['Data/' DataSetName '_data.mat']);

%reduce data size

data=reduce_data(data,N,rand_id);

%% Run Baseline

evalc('[result,data_mod] = apply_baseline_new(data);')
disp('Done Baseline');


%% Our Method %%
%initialize
[var] = initialize_admm_ndset(result,data_mod);
%Run IRLS based optimization
tic; [result1] = iterative_R3_irls_ndset(var,data_mod,MAXIT); tt=toc;
fprintf('Time takes for Optimization : %.2f sec\n',tt);
%fix rotation from Baseline
result1.R_out=result.R_out;
%Decompose Essentials into rotation and translation (apply baseline)
evalc('[result_admm_bs,data_admm] = apply_baseline_toadmm_noR_noT(data_mod,result1);')
disp('Done baseline on ADMM');


%% Result
disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
disp('Baseline');
 fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,100*result.errorF,100*result.errorF_median,100*result.Tran_err_mean,100*result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);
disp('ADMM');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,100*result_admm_bs.errorF,100*result_admm_bs.errorF_median,100*result_admm_bs.Tran_err_mean,100*result_admm_bs.Tran_err_median,result_admm_bs.Rot_err_mean,result_admm_bs.Rot_err_median);

