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

id = 1;
DataSetName=DataSetName1{id};

load(['Data/' DataSetName '_data.mat']);

%reduce data size

N = 50;
rand_id = 1;
data0=reduce_data(data,N,rand_id);

%% Run Baseline

[result,data_mod] = apply_baseline_new(data);

%%
evalc('[result,data_mod] = apply_baseline_new(data0);')
disp('Done Baseline');
disp('Baseline');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,data0.n,100*result.errorF,100*result.errorF_median,100*result.Tran_err_mean,100*result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);


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

result3=result;
result4=result1;
var0=var;

%%
N = 100;
rand_id = 1;
data0=reduce_data(data,N,rand_id);

[result,data_mod] = apply_baseline_new(data0);
[var] = initialize_admm_ndset(result,data_mod);

% var0=var;
AA = var.E_est./kron(var.lam,ones(3,3));
% AA = var.E_est;

AA(isnan(AA))=0;

cov = TMEnovel(AA,6,1);
% cov = TME(AA);

[U,S,V] = svd(cov);

% C = AA-U(1:6,:)'*U(1:6,:)*AA;
% C = U(7:end,:)'*U(7:end,:)*AA;
C = U(1:6,:)'*U(1:6,:)*AA;

scale = ((sum(C.*C,1)).^.5);

[~,Midx]=maxk(scale,20);

close all
plot(sort(scale))

rm_id = ceil(Midx/3);
rm_id = unique(rm_id);
data_rm = reduce_data_from_id(data0,rm_id);
rm_id

[result1,data_mod] = apply_baseline_new(data_rm);

fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,data0.n,100*result.errorF,100*result.errorF_median,100*result.Tran_err_mean,100*result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,data_rm.n,100*result1.errorF,100*result1.errorF_median,100*result1.Tran_err_mean,100*result1.Tran_err_median,result1.Rot_err_mean,result1.Rot_err_median);

%%
[result,data_mod] = apply_baseline_new(data_rm);
% evalc('[result,data_mod] = apply_baseline_new(data_rm);')
% disp('Done Baseline');
%
%
% %initialize
% [var] = initialize_admm_ndset(result,data_mod);
% %Run IRLS based optimization
% tic; [result1] = iterative_R3_irls_ndset(var,data_mod,MAXIT); tt=toc;
% fprintf('Time takes for Optimization : %.2f sec\n',tt);
% %fix rotation from Baseline
% result1.R_out=result.R_out;
% %Decompose Essentials into rotation and translation (apply baseline)
% evalc('[result_admm_bs,data_admm] = apply_baseline_toadmm_noR_noT(data_mod,result1);')
% disp('Done baseline on ADMM');
%


disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
disp('Baseline');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,100*result.errorF,100*result.errorF_median,100*result.Tran_err_mean,100*result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);
% disp('ADMM');
% fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,N,100*result_admm_bs.errorF,100*result_admm_bs.errorF_median,100*result_admm_bs.Tran_err_mean,100*result_admm_bs.Tran_err_median,result_admm_bs.Rot_err_mean,result_admm_bs.Rot_err_median);

%%

DataSetName=DataSetName1{id};
load(['Data/' DataSetName '_data.mat']);

N = 100;
rand_id = 1;
data=reduce_data(data,N,rand_id);

[result,data_mod] = apply_baseline_new(data);

%%
% [var0] = initialize_admm_ndset(result,data_mod);
AA = var0.E_est;
AA(isnan(AA))=0;
cov = TMEnovel(AA,6,4);
% cov = TME(AA);
[U,S,V] = svd(cov);
C = U(7:end,:)'*U(7:end,:)*AA;
scale = ((sum(C.*C,1)).^.5);
[~,Midx]=mink(scale,40);
% close all
% plot(sort(scale))
rm_id = ceil(Midx/3);
rm_id = unique(rm_id);
data_rm = reduce_data_from_id(data,rm_id);
rm_id

[result2,data_mod] = apply_baseline_new(data_rm);
disp('ID | rand_id | N | Mean Essential Error | Median Essential Error |  Mean Location Error | Median Location Error | Mean Rotation Error | Median Rotation Error |gm');
disp('Baseline');
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,data.n,100*result.errorF,100*result.errorF_median,100*result.Tran_err_mean,100*result.Tran_err_median,result.Rot_err_mean,result.Rot_err_median);
fprintf('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f\n',id,rand_id,data_rm.n,100*result2.errorF,100*result2.errorF_median,100*result2.Tran_err_mean,100*result2.Tran_err_median,result2.Rot_err_mean,result2.Rot_err_median);

%%
% close all
% plot(sort(scale))
[var1] = initialize_admm_ndset(result2,data_mod);
BB = var1.E_est;