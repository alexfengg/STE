function meas=run_measure_error(data)

meas.E_err=nan(data.n); meas.no_in=nan(data.n); meas.esti_mean=nan(data.n); meas.esti_median=nan(data.n);
meas.geo_mean=nan(data.n); meas.geo_median=nan(data.n);
for i=1:data.n
    for j=i+1:data.n
        if (i==j)||isempty(data.E_est{i,j})
            continue;
        end
        
        K1=data.K(:,:,i); K2=data.K(:,:,j);
        Egt=s3(data.E,i,j); E=data.E_est{i,j};
        pt=data.corr{i,j}; pt1=pt(:,:,1); pt2=pt(:,:,2);
        
        
       [E_err,no_in,esti_mean,esti_median,geo_mean,geo_median] = measure_error(E,Egt,pt1,pt2,K1,K2);
       if isnan(E_err)
           E_err=2;
       end
       
       meas.E_err(i,j)=E_err; 
       meas.no_in(i,j)=no_in;
       meas.esti_mean(i,j)=esti_mean;
       meas.esti_median(i,j)=esti_median;
       meas.geo_mean(i,j)=geo_mean;
       meas.geo_median(i,j)=geo_median;
       
    
        
    end
end

meas.E_err_vec=meas.E_err(:); meas.E_err_vec(isnan(meas.E_err_vec))=[];
meas.no_in_vec=meas.no_in(:); meas.no_in_vec(isnan(meas.no_in_vec))=[];
meas.esti_mean_vec=meas.esti_mean(:); meas.esti_mean_vec(isnan(meas.esti_mean_vec))=[];
meas.esti_median_vec=meas.esti_median(:); meas.esti_median_vec(isnan(meas.esti_median_vec))=[];
meas.geo_mean_vec=meas.geo_mean(:); meas.geo_mean_vec(isnan(meas.geo_mean_vec))=[];
meas.geo_median_vec=meas.geo_median(:); meas.geo_median_vec(isnan(meas.geo_median_vec))=[];
end