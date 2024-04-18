clc
clear
close all

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


%%
MEAN = [];
MEDIAN = [];
RES = [];
for i = 1
    dataname = data_names{i};
    save_file = ['result/RSR/',dataname,'_RSR.mat'];
    load(save_file)

    MEAN(i,:) = mean(R_err);
    MEDIAN(i,:) = median(R_err);

end

MEDIAN
MEAN