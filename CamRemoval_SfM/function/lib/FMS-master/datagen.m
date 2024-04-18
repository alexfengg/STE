function [data,SBSP] =  datagen(n,k,nin,nout,CV)

% generate inlier-outlier data
% Author: John Goes, 2014



%%%%%%%%% outliers
CV_sqrt = sqrtm(CV);
outs=CV_sqrt*random('normal',0,1,n,nout);


%%%%%%%%% inliers

    SBSP = random('normal',0,1,n,k);
    SBSP = orth(SBSP);
    SBSP2=SBSP*SBSP'*1/k;
    ins = SBSP2 * random('normal',0,1,n,nin);


% varin=var(ins')
% varout=var(outs')

data=[outs ins]';
data=data(randperm(size(data,1)),:);