% clc
clear
close all

% apply RSR (STE,TME,FMS,SFMS) methods to removing outlying camers for 
% photo tourism datasets. 20% of outlying columns of E are chosen as the
% outliers
% results are saved in "result/LUD+RSR/", the variable "result_RSR" contains
% the mean (median) errors of rotation and translation.

% add path
setpaths

%% 
data_names{1} = 'Alamo';
data_names{2} = 'Ellis_Island';
data_names{3} = 'Madrid_Metropolis';
data_names{4} = 'Montreal_Notre_Dame';
data_names{5} = 'Notre_Dame';
data_names{6} = 'NYC_Library';
data_names{7} = 'Piazza_del_Popolo';
data_names{8} = 'Piccadilly';
data_names{9} = 'Roman_Forum';
data_names{10} = 'Tower_of_London';
data_names{11} = 'Union_Square';
data_names{12} = 'Vienna_Cathedral';
data_names{13} = 'Yorkminster';
data_names{14} = 'Gendarmenmarkt';

mkdir('result/LUD+STE')
mkdir('result/LUD+STE+MC')
mkdir('result/LUD+TME')
mkdir('result/LUD+TME+MC')
mkdir('result/LUD+FMS')
mkdir('result/LUD+FMS+MC')
mkdir('result/LUD+SFMS')
mkdir('result/LUD+SFMS+MC')

opt.svtopt = 'normal';
opt.scaleopt = 'normal';

%%
for id=2
    %% load data
    DataSetName=data_names{id};
    load(['data/SfM_data/' DataSetName '_data.mat']);
    load(['data/bigE/',DataSetName,'_bigE.mat']);
    fprintf('Processing %s, n: %d\n', DataSetName, data.n)
    %% STE with Matrix Completion
    cov = STE(E2,6,3);
    [U,S,V] = svd(cov);
    C = U(:,1:6)*U(:,1:6)'*E2; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_STE,data_mod1,RemComp1,time_LUD_STE] = LUDpipline(data1);
    save(['result/LUD+STE+MC/',DataSetName,'_LUD_STE.mat'],'result_LUD_STE','time_LUD_STE')
    %% STE without Matrix Completion
    cov = STE(E1,6,3);
    [U,S,V] = svd(cov);
    C = U(:,1:6)*U(:,1:6)'*E1; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_STE,data_mod1,RemComp1,time_LUD_STE] = LUDpipline(data1);
    save(['result/LUD+STE/',DataSetName,'_LUD_STE.mat'],'result_LUD_STE','time_LUD_STE')
    %% TME with Matrix Completion
    cov = TME(E2);
    [U,S,V] = svd(cov);
    C = U(:,1:6)*U(:,1:6)'*E2; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_TME,data_mod_TME,RemComp_TME,time_LUD_TME] = LUDpipline(data1);
    save(['result/LUD+TME+MC/',DataSetName,'_LUD_TME.mat'],'result_LUD_TME','time_LUD_TME')
    %% TME without Matrix Completion
    cov = TME(E1);
    [U,S,V] = svd(cov);
    C = U(:,1:6)*U(:,1:6)'*E1; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_TME,data_mod_TME,RemComp_TME,time_LUD_TME] = LUDpipline(data1);
    save(['result/LUD+TME/',DataSetName,'_LUD_TME.mat'],'result_LUD_TME','time_LUD_TME')
    %% FMS with Matrix Completion
    U = fms_nosp(E2',6,opt);
    C = U(:,1:6)*U(:,1:6)'*E2; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_FMS,data_mod_FMS,RemComp_FMS,time_LUD_FMS] = LUDpipline(data1);
    save(['result/LUD+FMS+MC/',DataSetName,'_LUD_FMS.mat'],'result_LUD_FMS','time_LUD_FMS')
    %% FMS without Matrix Completion
    U = fms_nosp(E1',6,opt);
    C = U(:,1:6)*U(:,1:6)'*E1; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_FMS,data_mod_FMS,RemComp_FMS,time_LUD_FMS] = LUDpipline(data1);
    save(['result/LUD+FMS/',DataSetName,'_LUD_FMS.mat'],'result_LUD_FMS','time_LUD_FMS')
    %% SFMS with Matrix Completion
    E3 = E2./(repmat(sum(E2.*E2,2).^.5,1,size(E2,1)));
    E3(isnan(E3))=0; 
    U = fms_nosp(E3',6,opt);
    C = U(:,1:6)*U(:,1:6)'*E3; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_SFMS,data_mod_SFMS,RemComp_SFMS,time_LUD_SFMS] = LUDpipline(data1);
    save(['result/LUD+SFMS+MC/',DataSetName,'_LUD_SFMS.mat'],'cov','result_LUD_SFMS','time_LUD_SFMS')
    %% SFMS without Matrix Completion
    E3 = E1./(repmat(sum(E1.*E1,2).^.5,1,size(E1,1)));
    E3(isnan(E3))=0; 
    U = fms_nosp(E3',6,opt);
    C = U(:,1:6)*U(:,1:6)'*E3; % projection of data to the subspace
    Dist = ((sum(C.*C,1)).^.5); %  dist to the subspace
    [~,Midx]=maxk(Dist,round(data.n*0.2)); % choose outliers
    rm_id = ceil(Midx/3);
    rm_id = unique(rm_id);
    data1 = reduce_data_from_id(data,rm_id);
    [result_LUD_SFMS,data_mod_SFMS,RemComp_SFMS,time_LUD_SFMS] = LUDpipline(data1);
    save(['result/LUD+SFMS/',DataSetName,'_LUD_SFMS.mat'],'cov','result_LUD_SFMS','time_LUD_SFMS')

end

% delete path
rmpaths