import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import pickle
import random
from tqdm import tqdm
import numpy as np
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from IPython.display import clear_output, Image, display
import os
from scipy import sparse, io

plt.ion()
plt.show()

printFunctionNames = True


    
def loadData(inputDataset):
    """
    Load input dataset
    """
    if inputDataset == 'brainCIDR':
        path = '../input/brainCIDR/'
        df = pd.read_csv(f"{path}brainTags.csv", index_col = 0).T
        truth = pd.read_pickle(f'{path}truth.pkl')
    
    if inputDataset == 'pancreaticIsletCIDR':
        path = '../input/pancreaticIsletCIDR/'
        df = pd.read_csv(f"{path}pancreaticIsletTags.csv", index_col = 0).T
        truth = pd.read_pickle(f'{path}truth.pkl')
    
    if inputDataset == 'deng':
        path = '../input/deng/'
        df = pd.read_csv(f"{path}deng.csv", index_col = 0).T
        truth = pd.read_pickle(f'{path}truth.pkl')
        
    if inputDataset == 'celegans':
        path = '../input/celengans/'
        data = sparse.load_npz(f"{path}sparse_data.npz")
        data = data.todense()
        df1 = pd.DataFrame(data = data)
        df1.set_index(np.load(f"{path}cells.npy"), inplace=True)
        df1.columns = np.load(f"{path}genes.npy")
        return df1, None
    
    if inputDataset in [ 'sce10x_qc', 'sce2_qc', 'sce8_qc']:
        path = '../input/cellBench/'
        data = sparse.load_npz(f"{path}{inputDataset}.npz")
        data = data.todense()
        df = pd.DataFrame(data = data)
        df.set_index(np.load(f"{path}{inputDataset}_cells.npy"), inplace=True)
        df.columns = np.load(f"{path}{inputDataset}_genes.npy")
        truth = pd.read_pickle(f'{path}{inputDataset}_truth.pkl')
        
    return df, truth


