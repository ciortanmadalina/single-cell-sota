import sys
sys.path.append("..") # this adds to path parent directory in order to import utils file
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import random
from tqdm import tqdm
import numpy as np
from sklearn import metrics
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.mixture import GaussianMixture
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.decomposition import PCA
from sklearn.base import BaseEstimator

import gmm_utils
import utils

class LouvainModel(BaseEstimator):
    def __init__(self,  
                 minCellsPerGene=0, 
                 minGeneDispersion=0, 
                 log=True, 
                 scaler='standardScaleCells', 
                 pca_comp=10, 
                 doUmap=False, 
                 umap_comp=3, 
                 nb_neighbors=10, 
                 dataset = 'defaultDataset'):
        # Store copy of params optimisation
        self.params = locals()
        del self.params['self']

        self.minCellsPerGene = minCellsPerGene
        self.minGeneDispersion = minGeneDispersion
        self.log = log
        self.scaler = scaler
        self.pca_comp = pca_comp
        self.doUmap = doUmap
        self.umap_comp = umap_comp
        self.nb_neighbors = nb_neighbors
        self.dataset = dataset

    def preprocess(self,df):
        discreteDf = np.zeros(df.shape)
        discreteDf[np.where(df>0)] = 1
        genesToKeep = np.where(discreteDf.sum(axis = 0)>=self.minCellsPerGene )[0]
        df= df[df.columns[genesToKeep]]
        del discreteDf 

        # Remove genes which have a very low variance as they are expressed equally in all cells
        logDf = np.log1p(df)
        nonZeroMean = logDf.mean(axis = 0)
        nonZeroMean[nonZeroMean==0] = 1e-10
        dispersion =logDf.var(axis = 0)/nonZeroMean
        genesToKeep = np.where(dispersion>=self.minGeneDispersion )[0]
        df= df[df.columns[genesToKeep]]
        del logDf, nonZeroMean, dispersion

        if self.log:
            df = np.log1p(df)

        # scaling
        if self.scaler == 'none':
            scaledDf = df.values
        if self.scaler == 'standardScaleGenes':
            scaledDf = StandardScaler().fit_transform(df)
        if self.scaler == 'standardScaleCells':
            scaledDf = StandardScaler().fit_transform(df.T).T
        if self.scaler == 'robustScaleGenes':
            scaledDf = RobustScaler().fit_transform(df)
        if self.scaler == 'robustScaleCells':
            scaledDf = RobustScaler().fit_transform(df.T).T
        return scaledDf
    def fit(self, X, y):
        X = self.preprocess(X)
        self.pca= PCA(self.pca_comp)
        self.pca.fit(X)
        pass
    def predict(self, X):
        self.X = X.copy()
        X = self.preprocess(X)
        data = self.pca.transform(X)
        if self.doUmap:
            data = gmm_utils.getUmap(data, ncomp = self.umap_comp)
        clusters = gmm_utils.cluster_knn_louvain(data, neighbors = self.nb_neighbors)
        self.clusters = clusters
        return clusters
    
    def fitPredict(self, X):
        self.fit(X, None)
        return self.predict(X)
    
    def evaluate(self,y):
        """
        This method can be invoked after a predict has been performed
        as X and clusters are expected to be already computed
        """
        ev = gmm_utils.externalValidation(y, self.clusters)
        iv = gmm_utils.internalValidation(self.X, self.clusters)
        params = {**self.params, **ev, **iv}
        return params
    
    def scorer(self,clf, X, y_true):
        """
        This is the scored to be used by cross_val_score
        or grid searach CV
        """
        print(f'scorer {X.shape}, {y_true.shape}')
        from sklearn.metrics.cluster import adjusted_rand_score
        y_pred = clf.predict(X)
        return adjusted_rand_score(y_true, y_pred)
    
    def runAll(self):
        """
        This method will be used by the bayesian parameter search
        It returns input parameters enriched with validation scores
        """
        df, truth = utils.loadData(self.dataset)
        y = truth.clusters
        self.fitPredict(df)
        return self.evaluate(y)