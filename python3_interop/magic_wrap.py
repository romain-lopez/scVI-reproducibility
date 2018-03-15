import sys
import magic
import os
import numpy as np
import pandas as pd
import scipy.sparse



data_name = sys.argv[1]
n_cells = 0
if len(sys.argv) > 2:
    n_cells = int(sys.argv[2])

#load data
if data_name[-3:] == "npy":   
    X = np.load(data_name)
elif data_name[-3:] == "npz":
    X = scipy.sparse.load_npz(data_name).A

if n_cells == 0:
    n_cells = X.shape[0]
    
print("running MAGIC on " + str(n_cells) + " cells")
data = pd.DataFrame(X[:n_cells])
data.columns = [str(gene) for gene in data.columns.values]

# Load single-cell RNA-seq data
scdata = magic.mg.SCData(data, "sc-seq")
print(scdata)

scdata.run_magic(n_pca_components=20, random_pca=True, t=6, k=30, ka=10, epsilon=1, rescale_percent=99)

if len(sys.argv) == 2:
    np.save(data_name[:-4] + "_MAGIC.npy", scdata.magic.data.as_matrix())