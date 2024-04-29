import numpy as np
import cv2
import pydegensac
import multiprocess
import h5py
from tqdm import tqdm

def access_data(filePath):
    file = h5py.File(filePath)
    DATA = file['data']
    Keep = DATA['keep'][:, :]
    rows, cols = np.where(Keep > 0)
    edges = [[rows[i], cols[i]] for i in range(len(rows)) if rows[i] > cols[i]]
    data = []
    for edge in edges:
        i, j = edge
        ref = DATA['corr'][i, j]
        xi = DATA[ref][:, :, 0].T
        xj = DATA[ref][:, :, 1].T
        data.append([[i,j],xi,xj])
    return data

def compute_fundamental_matrices(data):
    edge, xi, xj = data
    F1, mask1 = cv2.findFundamentalMat(xi, xj, cv2.RANSAC, 0.5, 0.9999)
    F2, mask2 = pydegensac.findFundamentalMatrix(xi, xj, enable_degeneracy_check=True)
    F3, mask3 = pydegensac.findFundamentalMatrix(xi, xj, enable_degeneracy_check=False)
    return edge, F1, F2, F3, mask1, mask2, mask3

def parallel_pipeline(dataset_name):
    filePath = '../../CamRemoval_SfM/data/SfM_data/' + dataset_name + '_data.mat' 
    data = access_data(filePath)
    with multiprocess.Pool(processes=multiprocess.cpu_count()) as pool:
        results = list(tqdm(pool.imap(compute_fundamental_matrices, data), total=len(data)))
    return results

def store_hdf5(dataset_name,results):
    filename = f'result/{dataset_name}_ransac.h5'
    with h5py.File(filename, 'w') as f:
        for res in results:
            edge, F1, F2, F3, mask1, mask2, mask3 = res
            groupName = str(edge[1])+'-'+str(edge[0])
            group = f.create_group(groupName)
            group.create_dataset('F1',data=F1)
            group.create_dataset('F2',data=F2)
            group.create_dataset('F3',data=F3)
            group.create_dataset('mask1',data=mask1)
            group.create_dataset('mask2',data=mask2)
            group.create_dataset('mask3',data=mask3)

def main():
    DataSetNames = [
        'Alamo',
        'Ellis_Island',
        'Madrid_Metropolis',
        'Montreal_Notre_Dame',
        'Notre_Dame',
        'NYC_Library',
        'Piazza_del_Popolo',
        'Roman_Forum',
        'Tower_of_London',
        'Union_Square',
        'Vienna_Cathedral',
        'Yorkminster',
        'Gendarmenmarkt',
        'Piccadilly'
    ]
    for dataset_name in DataSetNames:
        print(f'Reading correspondence data: {dataset_name}')
        results = parallel_pipeline(dataset_name)
        # Process or store results here
        print('Done!')
        store_hdf5(dataset_name,results)

if __name__ == '__main__':
    main()
