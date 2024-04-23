clc
clear
close all

% apply LUD to photo tourism and save the results. also record the excution
% time for each run
% LUDpipline: original from "lib/SfM_CVPR2017_code/apply_baseline_new.m"
% results are saved in "result/LUD/", the variable "result_LUD" contains
% the mean (median) errors of rotation and translation.

% add path
setpaths

mkdir('result/LUD/')
%% Load data
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

for id=14
    DataSetName=data_names{id};

    load(['data/SfM_data/' DataSetName '_data.mat']);
    fprintf('Processing %s, n: %d\n', DataSetName, data.n)
    
    [result_LUD,data_mod_LUD,RemComp,time_LUD] = LUDpipline(data);
    time_LUD
    result_LUD
    
    save(['result/LUD/',data_names{id},'_LUD.mat'],'result_LUD','data_mod_LUD','RemComp','time_LUD');
    
end

rmpaths