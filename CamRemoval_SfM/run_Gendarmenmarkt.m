% clc
clear
close all

% apply RSR to remove the cameras of Gendarmenmarkt dataset with 45%
% outlying columns

% add path
setpaths

%% Load data

DataSetName='Gendarmenmarkt';

load(['data/SfM_data/',DataSetName,'_data.mat']);
load(['data/bigE/',DataSetName,'_bigE.mat'])

fprintf('Processing %s, n: %d\n', DataSetName, data.n)

opt.scaleopt = 'normal';
opt.svdopt = 'normal';
%% FMS, the removal rate is 45%
U = fms(E2',6,opt);
C = U(:,1:6)*U(:,1:6)'*E2; % projection of data to the subspace
Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
[~,Midx]=maxk(Dist,round(data.n*0.4)); % choose outliers
rm_id = unique(ceil(Midx/3));
data1 = reduce_data_from_id(data,rm_id);
[result_FMS,~,~,time_FMS] = LUDpipline(data1);

%% STE
cov = STE(E2,6,3);
[U,S,V] = svd(cov);
C = U(:,1:6)*U(:,1:6)'*E2; % projection of data to the subspace
Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
[~,Midx]=maxk(Dist,round(data.n*0.45)); % choose outliers
rm_id = unique(ceil(Midx/3));
data1 = reduce_data_from_id(data,rm_id);
[result_STE,~,~,time_STE] = LUDpipline(data1);

%% TME
cov = TME(E2);
[U,S,V] = svd(cov);
C = U(:,1:6)*U(:,1:6)'*E2; % projection of data to the subspace
Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
[~,Midx]=maxk(Dist,round(data.n*0.45)); % choose outliers
rm_id = unique(ceil(Midx/3));
data1 = reduce_data_from_id(data,rm_id);
[result_TME,~,~,time_TME] = LUDpipline(data1);

%% SFMS
F = E2./(repmat(sum(E2.*E2,2).^.5,1,size(E2,1)));
F(isnan(F))=0;
U = fms(F',6,opt);
C = U(:,1:6)*U(:,1:6)'*F; % projection of data to the subspace
Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
[~,Midx]=maxk(Dist,round(data.n*0.45)); % choose outliers
rm_id = unique(ceil(Midx/3));
data1 = reduce_data_from_id(data,rm_id);
[result_SFMS,~,~,time_SFMS] = LUDpipline(data1);

%%
save_dir = 'result/Gendarmenmarkt-45/';
if ~exist(save_dir)
    mkdir(save_dir)
end

save([save_dir,DataSetName,'_45.mat'],'result_TME','result_STE', ...
    'result_SFMS','time_TME','time_STE','time_SFMS')

time_STE
time_TME
time_FMS
time_SFMS

result_STE
result_TME
result_FMS
result_SFMS



