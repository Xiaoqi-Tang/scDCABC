import os
import pandas as pd
import numpy as np
import networkx as nx
import scanpy as sc
from sklearn.cluster import KMeans
import community.community_louvain as community_louvain

from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist, squareform
#import community
import leidenalg as la
import igraph as ig
import community.community_louvain as community_louvain
import dca_main
import gc
import time

def rm1_g(x1, x2, a):
    class_indices = np.where(x2 == a)[0]
    #ma = x1[:, class_indices] #numpy数组的索引方法
    ma = x1.iloc[:, class_indices] #数据框的索引方法
    return ma

def rm1(x1, x2, a):
    class_indices = np.where(x2 == a)[0]
    ma = x1.iloc[class_indices, :]
    return ma

#def list_frame(x):
#    data_frame = pd.DataFrame(columns=['v1', 'v2'])
#    path_counter = 0
#    for path in x.keys():
#        path_counter += 1
#        temp_data = x[path]
#        temp_data['v2'] += path_counter - 1
#        data_frame = pd.concat([data_frame, temp_data], ignore_index=True)
#   return data_frame

def list_frame(x):
    data_frame = pd.DataFrame(columns=['v1', 'v2'])
    path_counter = 0
    max_length = 0  
    longest_path = ""  

    for path in x.keys():
        path_counter += 1
        temp_data = x[path]
        temp_data['v2'] += path_counter - 1
        data_frame = pd.concat([data_frame, temp_data], ignore_index=True)
        
        path_numbers = path.split('->')
        num_count = len(path_numbers)
        
        if num_count > max_length:
            max_length = num_count
            longest_path = path

    return data_frame

def gene_c(x1):
    np.random.seed(123)
    result = KMeans(n_clusters=5, n_init=20).fit(x1.T)
    clust1 = pd.DataFrame(result.labels_, columns=['cluster'])
    unique_clust = np.sort(clust1['cluster'].unique())

    min_mean_dist = np.inf
    matrix_g = None

    for a in unique_clust:
        rm11_g = rm1_g(x1, clust1['cluster'].values, a)
        #if rm11_g.shape[0] < 2 or np.all(rm11_g.sum(axis=1) == 0) or rm11_g.shape[1] < 2 or not isinstance(rm11_g, pd.DataFrame):
        if len(rm11_g) < 3 or np.all(rm11_g.sum(axis=1) == 0) or rm11_g.shape[1] < 3 or not isinstance(rm11_g, pd.DataFrame):
            continue
        else:
            new = cell_c(rm11_g)
            clust2 = new.iloc[:, 1]
            unique_clust2 = np.sort(clust2.unique())

            mean_ddd = []
            for b in unique_clust2:
                rm11_cc = rm1(rm11_g, clust2.values, b)
                dist_matrix = pairwise_distances(rm11_cc, metric='euclidean')
                mean_d = dist_matrix.mean()
                mean_ddd.append(mean_d)
                mean_dist = np.nanmean(mean_ddd)
            if mean_dist < min_mean_dist:
                min_mean_dist = mean_dist
                matrix_g = rm11_g

    return matrix_g

# Walktrap
def cell_c(x):
    latent3 = dca(x.T) 
    latent = pd.DataFrame(latent3) 
    np.random.seed(123)
    adata = sc.AnnData(X=latent)
    adata = sc.AnnData(X=x)
    adata.layers["logcounts"] = adata.X.copy()
    sc.pp.neighbors(adata, n_neighbors=min(adata.n_obs - 1, 20), use_rep="X")

    adj_matrix = adata.obsp["connectivities"]
    g = ig.Graph.Weighted_Adjacency(adj_matrix, mode='UNDIRECTED', attr='weight')


    walktrap = g.community_walktrap(weights='weight').as_clustering()
    clust = walktrap.membership
    clust1 = pd.DataFrame({'v1': x.index, 'v2': clust})
    return clust1

def re(x1, a, data=None, path=[]):
    if data is None:
        data = {}
    clust1 = cell_c(gene_c(x1))
    k = clust1['v2'].max()

    for b in range(0, k+1):
        new_path = path + [b]
        if k == 0:
            data["->".join(map(str, [a] + path))] = clust1[clust1['v2'] == b]
        elif clust1[clust1['v2'] == b].shape[0] < 6:
            data["->".join(map(str, [a] + new_path))] = clust1[clust1['v2'] == b]
        else:
            x11 = rm1(x1, clust1['v2'].values, b)
            data = re(x11, a, data, new_path)
    return data

def main(x):
    clust = cell_c(x)
    unique_clust = np.sort(clust['v2'].unique())
    data_list = {}
    for a in unique_clust:
        rm11_c = rm1(x, clust['v2'].values, a)
        data_list.update(re(rm11_c, a))
    data = list_frame(data_list)
    return data

def dca(x):
    print(x.shape)
    args = dca_main.parse_args()
    args.input = x
    args.outputdir = "C:/Users/W10/Desktop/latent"
    args.batchsize = 16  
    try:
        import tensorflow as tf
    except ImportError:
        raise ImportError('DCA requires TensorFlow v2+. Please follow instructions'
                          ' at https://www.tensorflow.org/install/ to install'
                          ' it.')
    # import tf and the rest after parse_args() to make argparse help faster
    import train
    #train.train_with_args(args)
    latent2 = train.train_with_args(args)
    gc.collect()
    return latent2



start_time = time.time()
matrix = pd.read_csv("D:/data/GSE36522/GSE36522.csv", index_col=0, header=0) #列是细胞 行是基因
matrix = matrix.apply(pd.to_numeric)
matrix1 = matrix.T 
result = main(matrix1)
print(result)
label = pd.read_csv("C:/Users/W10/Desktop/DD/GSE36522/label.csv", dtype=str)
unique_values = result['v2'].unique()

data = pd.DataFrame()
for i in unique_values:
    sub_result = result[result['v2'] == i].copy()
    sub_result['v2'] = np.nan
    for j in range(len(sub_result)):
        element = sub_result['v1'].iloc[j]
        match_idx = label[label['V1'] == element].index
        if not match_idx.empty:
            label_value = label.loc[match_idx[0], 'V2']
            #sub_result['v2'].iloc[j] = label_value
            sub_result.loc[sub_result.index[j], 'v2'] = label_value  # 使用loc方法进行赋值
    label_count = sub_result['v2'].value_counts()
    max_label = label_count.idxmax()
    sub_result['v2'] = max_label
    data = pd.concat([data, sub_result], ignore_index=True)

data_sorted = data.set_index('v1').reindex(label['V1']).reset_index()
data_sorted.to_csv("C:/Users/W10/Desktop/scDCABC.csv", index=False)


end_time = time.time()

run_time_minutes = (end_time - start_time) / 60
print(f"time: {run_time_minutes:.2f} ")