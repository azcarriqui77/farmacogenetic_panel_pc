import numpy as np
import pandas as pd
import random
import time
import os
import multiprocessing as mp

def dataLoading(path, features_in_columns):
    data = pd.read_csv(path)
    if features_in_columns == False:
        data = data.transpose()
    target = data["target"].astype('category').tolist()
    df = data.drop(columns=["target"]).sample(frac=1, axis=1)
    return df, target
    
def splitDataByColumn(df, n_clusters):
    # Select the lowest value for number of processors to use
    clusters = min(mp.cpu_count(), n_clusters, df.shape[1])
    
    # Split dataframe by columns and append each one to a list
    columns = df.columns.tolist()
    columns_split = np.array_split(columns, clusters)
    dfs_particionados = []
    for subset in columns_split:
        dfs_particionados.append(df[list(subset)])
    
    return dfs_particionados, clusters
            
