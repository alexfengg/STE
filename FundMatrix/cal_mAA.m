clc
clear
close all

% calculate mAA(10) from the saved results

data_names = {
    'Alamo';
    'Ellis_Island';
    'Madrid_Metropolis';
    'Montreal_Notre_Dame';
    'NYC_Library';
    'Notre_Dame';
    'Piazza_del_Popolo';
    'Piccadilly';
    'Roman_Forum';
    'Tower_of_London';
    'Union_Square';
    'Vienna_Cathedral';
    'Yorkminster';
    'Gendarmenmarkt'
    };

df = 0.01;
thr = [1:10];

mAA1 = zeros(14,4);
mAA = cell(1,14);

for m=2
    dataName = data_names{m};
    
    save_file = ['result/RSR/',dataName,'_RSR.mat'];
    load(save_file)
    R=R_err;
    N = size(R,1);
    for i = 1:4
        tmp = zeros(10,1);
        for j=1:size(thr,2)
            tmp(j) = sum(R(:,i)<thr(j))/N;
        end
        mAA{m}(i) = mean(tmp);

    end
    mAA1(m,:) = mAA{m};
end

mAA1