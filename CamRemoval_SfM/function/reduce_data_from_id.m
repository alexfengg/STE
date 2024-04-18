function [data,id] = reduce_data_from_id(data,id_rm)

data.n = data.n-length(id_rm);

id=id_rm;
id3=[3*id-2;3*id-1;3*id]; id3=id3(:);

if size(data.t,2)==1
    data.t(:,:,id) = [];
else
    data.t(:,id) = [];
end

data.R(:,:,id) = [];
data.keep(:,id) = [];
data.keep(id,:) = [];

data.E(:,id3) = [];
data.E(id3,:) = [];

data.E_gt(:,id3) = [];
data.E_gt(id3,:) = [];

data.G_gt(:,id3) = [];
data.G_gt(id3,:) = [];

data.Hmat(:,id3) = [];
data.Hmat(id3,:) = [];

data.E_est(:,id) = [];
data.E_est(id,:) = [];

data.corr(:,id) = [];
data.corr(id,:) = [];

data.tijGT(:,id) = [];
data.tijGT(id,:) = [];

data.AdjMat(:,id) = [];
data.AdjMat(id,:) = [];

data.K(:,:,id) = [];
data.Focal_gt(id) = [];
data.imgW(id) = [];
data.imgL(id) = [];
data.W_x_full(:,id) = []; 
data.W_y_full(:,id) = [];
data.W_mask_full(:,id) = [];
data.CompInds(id) = [];


end