function [data] = reduce_data_from_id_edges(data,id_rm)

for i=id_rm
    for j=id_rm
        data.AdjMat(i,j)=0;
        data.Hmat([3*i-2,3*i-1,3*i],[3*j-2,3*j-1,3*j]) = 0;
        
        data.E_est{i,j} = {};
        data.keep(i,j)=0;
    end
end

end