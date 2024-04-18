function [S_pruned,R,t,K,W_x_pruned,W_y_pruned,W_mask_pruned,cams,pts3D]=create_3D_initial(data,result,norm_lim,max_th,av_th,ext)
%clc
addpath(genpath('../../ToSoumyadip'));

 %data=reduce_data(data,100,1);
% 
% [result,data]=apply_baseline_new(data);

%% convert to correct format
% J=diag([1 1 -1]); 
% for i=1:data.n
%     R(:,:,i)=J*data.R(:,:,i)*J;
%     t(:,:,i)=-J*data.R(:,:,i)*data.t(:,i);
%     K(:,:,i)=eye(3); K(2,2,i)=data.K(2,2,i); K(1,1,i)=data.K(1,1,i);
% end

J=diag([1 1 -1]); 
for i=1:data.n
    R(:,:,i)=J*result.R_out{i}*J;
    t(:,:,i)=-J*result.R_out{i}*result.T_out(:,i);
    K(:,:,i)=eye(3); K(2,2,i)=-data.K(2,2,i); K(1,1,i)=data.K(1,1,i);
end


%% create imagelist

for kk = 1:data.n
    ImageList(1,kk).FName = num2str(kk);
    ImageList(1,kk).f = data.K(1,1,kk);
    ImageList(1,kk).Width = 2*data.K(1,3,kk);
    ImageList(1,kk).Height = 2*data.K(2,3,kk);
    ImgWidth(1,kk)=2*data.K(1,3,kk);
    ImgLength(1,kk)=2*data.K(2,3,kk);
end

%% create correspondence
CorrPts=data.corr;
CoorPtsCntred = cell(data.n,data.n);
CorrPtsNums=zeros(data.n);

for kk = 1:data.n
    for ss = kk+1:data.n
        if ~isempty(CorrPts{kk,ss})
            CoorPtsCntred{kk,ss} = CorrPts{kk,ss};
            CoorPtsCntred{kk,ss}(1,:,1) = CorrPts{kk,ss}(1,:,1)-ImgWidth(kk)/2;
            CoorPtsCntred{kk,ss}(2,:,1) = -(CorrPts{kk,ss}(2,:,1)-ImgLength(kk)/2);
            CoorPtsCntred{kk,ss}(1,:,2) = CorrPts{kk,ss}(1,:,2)-ImgWidth(ss)/2;
            CoorPtsCntred{kk,ss}(2,:,2) = -(CorrPts{kk,ss}(2,:,2)-ImgLength(ss)/2);
            CorrPtsNums(kk,ss)=size(CoorPtsCntred{kk,ss},2);
        end
    end
end

AllPairsCompact.pts = CoorPtsCntred;
AllPairsCompact.pts_num = CorrPtsNums;

%% initalize 3D
W_x_ff=data.W_x_full; W_y_ff=data.W_y_full; W_mask_ff=data.W_mask_full;
%track=create_track(data); W_x_ff=track(:,:,1); W_y_ff=track(:,:,2); W_mask_ff=(track(:,:,2)~=0);
ImgSizes = [ImgWidth; ImgLength];


prm.nonlin = true;
S_init = TriangChains(W_x_ff, W_y_ff, W_mask_ff, R, t, K, ImgSizes, prm);

%% remove nan structure points
ff_valid = ~any(isnan(S_init),1);
S = S_init(:,ff_valid);
W_x = W_x_ff(ff_valid, :);
W_y = W_y_ff(ff_valid, :);
W_mask = W_mask_ff(ff_valid, :);

%% Remove large norm points
S_origin = median(S,2);
S_ctrd = bsxfun(@minus, S, S_origin);

SnormsVec = norms(S_ctrd,2,1);
hist(SnormsVec(SnormsVec<100),100);
normLim = norm_lim;%quantile(SnormsVec,0.98); %norm_lim
S = S(:,SnormsVec<normLim);
W_x = W_x(SnormsVec<normLim, :);
W_y = W_y(SnormsVec<normLim, :);
W_mask = W_mask(SnormsVec<normLim, :);
%% display structure
S_origin = median(S,2);
S_ctrd = bsxfun(@minus, S, S_origin);


%% calculate reprojection errors
[W_err, mse_err, mean_err, num_reproj] = ComputeReprojectionErr_Fast(S, R, t, K, W_x, W_y, W_mask, 1);
W_err_avg = sum(W_err,2)./sum(W_mask,2);
%% select invalid/high reprojection error chains
reproj_thresh_max = quantile(W_err(:),max_th);%max_th=0.995;
reproj_thresh_avg = quantile(W_err_avg(:),av_th);%av_th=0.95;
min_chain_length = 2;


ff_chains_bad = any(isnan(W_err),2) | (max(W_err,[],2)>reproj_thresh_max) | (W_err_avg>reproj_thresh_avg) | (sum(W_mask,2)<min_chain_length);

fprintf('Removing %d out of %d chains, %d chains remain...\n', nnz(ff_chains_bad), length(ff_chains_bad), nnz(~ff_chains_bad));

W_mask_pruned = W_mask(~ff_chains_bad,:);
sum(sum(W_mask_pruned,2)>2)

W_x_pruned = W_x(~ff_chains_bad,:);
W_y_pruned = W_y(~ff_chains_bad,:);
S_pruned = S(:,~ff_chains_bad);

S_origin_pruned = median(S_pruned,2);
S_ctrd_pruned = bsxfun(@minus, S_pruned, S_origin_pruned);

figure;
plot3(S_ctrd_pruned(1,:), S_ctrd_pruned(2,:), S_ctrd_pruned(3,:), 'b.', 'markersize',1);
axis([-1 1 -1 1 -1 1]*normLim);


%% alternate camera file
cam_id=fopen(['../../sba-1.6/demo/cam1_real' ext '.txt'],'w+');
for k=1:data.n
    quat = rotm2quat(R(:,:,k));
    t_tr=t(:,1,k);%-R(:,:,k)'*t(:,1,k);
    fprintf(cam_id,'%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ... 
        K(1,1,k),K(1,3,k),K(2,3,k),K(2,2,k)/K(1,1,k),0,quat(1),quat(2),quat(3),quat(4),t_tr(1),t_tr(2),t_tr(3));
end
fclose(cam_id);

% % cal_id=fopen('/fs/jacobsface/SfM/sba-1.6/demo/calib_real.txt','w+');
% % for k=1:data.n
% %     fprintf(cal_id,'%.6f %.6f %.6f %.6f %.6f\n',K(1,1,k),K(1,2,k),K(2,2,k),K(1,3,k),K(2,3,k));
% % end
% % fclose(cal_id);

%% write 3d points file
X3d_id=fopen(['../../sba-1.6/demo/X3d_real' ext '.txt'],'w+');
for i=1:size(W_mask_pruned,1)
    
    X=S_pruned(:,i);
    tr=cat(3,full(W_x_pruned(i,:)),full(W_y_pruned(i,:)));
    frame=find(W_mask_pruned(i,:)~=0);
    co_od=squeeze(tr(:,frame,:));
    
    fprintf(X3d_id,'%.6f %.6f %.6f %d ',X(1),X(2),X(3),size(co_od,1));
    for tt=1:size(co_od,1)
        fprintf(X3d_id,'%d %.6f %.6f ',frame(tt)-1,co_od(tt,1),co_od(tt,2));
    end
    fprintf(X3d_id,'\n');   
end
fclose(X3d_id);

%% write unix code part to do BA
out_file=['../../sba-1.6/demo/result' ext '.txt'];
delete(out_file);
cur_dir=pwd;
cd ../../sba-1.6/demo/
unix(['./eucsbademo cam1_real' ext '.txt X3d_real' ext '.txt >> result' ext '.txt']);
cd(cur_dir);

%% colect output


o_id=fopen(out_file,'r');

for i=1:7
    fgetl(o_id);
end


for i=1:data.n
    par(i,:)=str2num(fgetl(o_id));
end
cams=par(:,6:end);

for i=1:3
    fgetl(o_id);
end

pt=textscan(o_id,'%f %f %f');
pts3D=[pt{1} pt{2} pt{3}];

fclose(o_id);

end


