  function [X3d]=create_input_BA(data)

%create track
track=create_track(data); %N_track*n_cam*2

%% delete some tracks
idx=sum(track(:,:,1)~=0,2);
id_del=(idx==2);
track(id_del,:,:)=[];


X3d=initialize_X3d(data,track); %3*N_track format

%% write camera file
% cam_id=fopen('cam.txt','a');
% for k=1:data.n
%     K=data.K(:,:,k);
%     fx=K(1,1); xo=K(1,3); yo=K(2,3); rat=K(2,2)/K(1,1); skew=0;
%     quat = rotm2quat(data.R(:,:,k));
%     fprintf(cam_id,'%.6f %.6f %.6f %.6f %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ... 
%         fx,xo,yo,rat,skew,quat(1),quat(2),quat(3),quat(4),data.t(1,1,k),data.t(2,1,k),data.t(3,1,k));
% end
% fclose(cam_id);

%% alternate camera file
cam_id=fopen('cam1_real.txt','w+');
for k=1:data.n
    quat = rotm2quat(data.R(:,:,k));
    fprintf(cam_id,'%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ... 
        quat(1),quat(2),quat(3),quat(4),data.t(1,1,k),data.t(2,1,k),data.t(3,1,k));
end
fclose(cam_id);

cal_id=fopen('calib_real.txt','w+');
for k=1:data.n
    fprintf(cal_id,'%.6f %.6f %.6f %.6f %.6f\n',data.K(1,1,k),data.K(1,2,k),-data.K(2,2,k),data.K(1,3,k),data.K(2,3,k));
end
fclose(cal_id);

%% write 3d points file
X3d_id=fopen('X3d_real.txt','w+');
for i=1:size(track,1)
    
    X=X3d(:,i);
    tr=track(i,:,:);
    frame=find(tr(:,:,1)~=0);
    co_od=squeeze(tr(:,frame,:));
    
    fprintf(X3d_id,'%.6f %.6f %.6f %d ',X(1),X(2),X(3),size(co_od,1));
    for t=1:size(co_od,1)
        fprintf(X3d_id,'%d %.6f %.6f ',frame(t)-1,co_od(t,1),co_od(t,2));
    end
    fprintf(X3d_id,'\n');   
end
fclose(X3d_id);