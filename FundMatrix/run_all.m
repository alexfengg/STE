clc
clear
close all

% running robust subspace methods (STE,TME,FMS,SFMS) for estimating
% fundamental matrix of the photo tourism datasets

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

setpaths
%%
for m = 2
    dataName = data_names{m};
    % load(['data/SfM_data/',dataName,'_data.mat'])
    load(['data/sfm_clean/',dataName,'_data_clean.mat'])

    N = size(datum,2);
    N = 100;

    message = ['Dataset: ',dataName,'. Total number of samples: ',num2str(N)];
    disp(message)
    R_err = parforComupte(datum,N,message);
    disp('Finished! The median errors:')
    median(R_err)

    save_loc = 'result/RSR/';
    if ~exist(save_loc,'dir')
        mkdir(save_loc)
    end
    save([save_loc,dataName,'_RSR.mat'],'R_err')
end

rmpaths
%%

function R_err = parforComupte(datum,N,message)

R_err = zeros(N,4);
D = parallel.pool.DataQueue;
h = waitbar(0, message);
afterEach(D, @nUpdateWaitbar);
p=1;

opt.svtopt = 'normal';
opt.scaleopt = 'normal';

parfor i = 1:N
    send(D, i);
    % fprintf('Dataset: %s (%d/14), %d / %d\n', dataName, m,i,N);
    df = datum{i};
    R0 = df.R0;
    x1 = df.points(1:2,:);
    x2 = df.points(3:4,:);
    K1 = df.K(:,1:3);
    K2 = df.K(:,4:6);
    [x1n,T1] = normalize_points(x1);
    [x2n,T2] = normalize_points(x2);
    X = repmat(x2n,3,1).*repelem(x1n,3,1);

    %% STE
    F1 = STE_FundMat(X);
    F1 = T1'*F1*T2;

    %% TME
    cov = TME(X);
    [U,S,V] = svd(cov);
    F = reshape(V(:,end),[3,3])';
    [U,S,V] = svd(F);
    S(3,3)=0;
    F2 = U*S*V';
    F2 = T1'*F2*T2;

    %% SFMS
    [L1,L2] = fms(X',8,opt);
    F = reshape(L2,[3,3])';
    [U,S,V] = svd(F);
    S(3,3)=0;
    F3 = U*S*V';
    F3 = T1'*F3*T2;

    %% FMS
    [L1,L2] = fms_nosp(X',8,opt);
    F = reshape(L2,[3,3])';
    [U,S,V] = svd(F);
    S(3,3)=0;
    F4 = U*S*V';
    F4 = T1'*F4*T2;

    %% SAVE RESULTS

    [R1,t1] = relpos2(F1,x1,x2,K1,K2);
    [R2,t2] = relpos2(F2,x1,x2,K1,K2);
    [R3,t3] = relpos2(F3,x1,x2,K1,K2);
    [R4,t4] = relpos2(F4,x1,x2,K1,K2);

    a1 = abs(acos((trace(R1*R0')-1)/2))/pi*180;
    a2 = abs(acos((trace(R2*R0')-1)/2))/pi*180;
    a3 = abs(acos((trace(R3*R0')-1)/2))/pi*180;
    a4 = abs(acos((trace(R4*R0')-1)/2))/pi*180;
    R_err(i,:) = [a1,a2,a3,a4];
end
close(h)

    function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end

end