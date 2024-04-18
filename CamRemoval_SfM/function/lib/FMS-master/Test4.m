%% Example script for running the FMS algorithm
%
% Author: Tyler Maunu, 2014

d=5;
D=100;
N_in=100;
N_out=100;
N=N_in+N_out;
noise=1e-8;
options=struct;

options.scaleopt='normal';      %normal or log
options.p=1;                    %0<p<=1
options.initopt='pca';          %random or pca
options.svdopt='randomized';    %normal or randomized
options.maxiter=100;            %default
options.epsilon=10^-10;         %default


[X,direction]=datagen(D,d,N_in,N_out,eye(D)/D);
X=X+random('normal',0,noise,N,D);


L=fms(X,d,options);
calc_sdist(L,direction)



