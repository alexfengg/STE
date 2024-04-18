function [data,id] = reduce_data(data,N,rand_id)
% rng(1);
rng('default');

gap=floor((data.n-N)/5);

st=1:gap:data.n; en=N:gap:data.n;
exp_no=min(numel(st),numel(en));
st=st(1:exp_no); en=en(1:exp_no);
rnid=randperm(exp_no); 
st=st(rnid); en=en(rnid);
 
id_ini=st(rand_id):en(rand_id); 
%id_ini=1:N;
keep=data.keep(id_ini,id_ini);
%% prune already disconnected edges
id_del=(sum(keep,1)==0);
id=id_ini(~id_del);


id3=[3*id-2;3*id-1;3*id]; id3=id3(:);

data.t=data.t(:,:,id);
data.R=data.R(:,:,id);
data.n=N;
data.keep=data.keep(id,id);
data.E=data.E(id3,id3);
data.E_gt=data.E_gt(id3,id3);
data.G_gt=data.G_gt(id3,id3);
data.Hmat=data.Hmat(id3,id3);
data.E_est=data.E_est(id,id);
try data.corr=data.corr(id,id);
catch
    data.corr=[];
end
try data.tijGT=data.tijGT(id,id);
catch
    data.tijGT=[];
end
data.AdjMat=data.AdjMat(id,id);
data.K=data.K(:,:,id);
data.Focal_gt=data.Focal_gt(id);
data.imgW=data.imgW(id);
data.imgL=data.imgL(id);

%process track
data.W_x_full=data.W_x_full(:,id); data.W_y_full=data.W_y_full(:,id);
data.W_mask_full=data.W_mask_full(:,id);

kp_idx=(sum(data.W_mask_full,2)>1);
data.W_x_full=data.W_x_full(kp_idx,:); data.W_y_full=data.W_y_full(kp_idx,:);
data.W_mask_full=data.W_mask_full(kp_idx,:);


data.n=size(data.t,3);
data.CompInds=data.CompInds(id);

end