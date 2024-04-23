% clc
clear
close all

% generate the big E matrix with/without matrix completion (MC)
% big E without MC: E1 contains (i,j)-th block Eij: scaled estimated, lam_ij*Eij
% Eij: estimated essential matrix fro image i,j
% lam_ij: scaler <Eij,Eij^LUD>/||Eij^LUD||^2_F, which is an initial value
% of the scaler used in [1]

% With MC: E1 contains many zero blocks, filling zeros blocks by running MC
% algorithm, Singular Value Thresholding (SVT), store the result as E2

% [1]: Soumyadip Sengupta, Tal Amir, Meirav Galun, Tom Goldstein, 
% David W. Jacobs, Amit Singer, Ronen Basri, "A New Rank Constraint on 
% Multi-view Fundamental Matrices, and its Application to Camera
% Location Recovery"

% [2]: Jian-Feng Cai, Emmanuel J. Candes, Zuowei Shen, "A Singular Value 
% Thresholding Algorithm for Matrix Completion"

% add path
setpaths

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

if ~exist('data/bigE','dir')
    mkdir('data/bigE')
end

for id=1:14
    %% load data
    DataSetName=data_names{id};
    load(['data/SfM_data/' DataSetName '_data.mat']);
    load(['result/LUD/',DataSetName,'_LUD.mat'])
    fprintf('Processing %s, n: %d\n', DataSetName, data.n)
    
    %% evaluate the scaler lam_ij
    lam = generate_initial_lambda(result_LUD,data_mod_LUD,data,RemComp);
    E0 = generateFstacked(data,lam);
    E0(isnan(E0))=0;
    E1 = E0;
    %% MC
    E1(isnan(E1))=0;
    Omega = E1~=0;
    
    % adjust the stepsize for some datasets smaller to avoid divergence
    if ismember(id,[8,9,10,11,14])
        delta = 0.1*size(E1,1)^2/sum(Omega(:)); % step size
    else
        delta = [];
    end
    X1 = SVT(E1,Omega,delta);
    X = X1;
    n = size(X,1)/3;
    omega1 = Omega(1:3:end,1:3:end);
    for i=1:n
        for j=1:n
            if omega1(i,j)==0
                TMP = X(3*i-2:3*i,3*i-2:3*i);
                [U,S,V]=svd(TMP);
                X(3*i-2:3*i,3*i-2:3*i) = U*diag([1,1,0])*U';
            end
        end
    end
    for i=1:n
        X(3*i-2:3*i,3*i-2:3*i) = 0;
    end

    E2 = E1+X.*(1-Omega);
    save(['data/bigE/',DataSetName,'_bigE.mat'],'E1','E2')

end
