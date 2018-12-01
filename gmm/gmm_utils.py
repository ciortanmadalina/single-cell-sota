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
from sklearn.datasets.samples_generator import make_blobs
import itertools
from scipy.spatial.distance import cdist
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.metrics.cluster import adjusted_rand_score
import umap
import os
plt.ion()
plt.show()

printFunctionNames = True

if printFunctionNames:
    print('elbowAnalysis')
def elbowAnalysis(X, numberOfClusters):
    distortions = []

    for k in tqdm(numberOfClusters):
        kmeanModel = KMeans(n_clusters=k)
        kmeanModel.fit(X)
        distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])

    plt.plot(numberOfClusters, distortions, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Distortion')
    plt.title('The Elbow Method showing the optimal k')
    plt.show()
    
    
if printFunctionNames:
    print('silhouetteAnalyis')
def silhouetteAnalyis (X, numberOfClusters):
    silhouette_score_values=[]
    for i in tqdm(numberOfClusters):
        classifier=KMeans(i,init='k-means++', n_init=10, max_iter=300, tol=0.0001, verbose=0, random_state=None, copy_x=True)
        classifier.fit(X)
        labels= classifier.predict(X)
        silhouette_score_values.append(metrics.silhouette_score(X,labels ,metric='euclidean', sample_size=None, random_state=None))

    plt.plot(numberOfClusters, silhouette_score_values)
    plt.title("Silhouette score values vs Numbers of Clusters ")
    plt.show()

    Optimal_NumberOf_Components=numberOfClusters[silhouette_score_values.index(max(silhouette_score_values))]
    print( "Optimal number of components is:", Optimal_NumberOf_Components)
    
    
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
    return df, truth


def optimalNbClustersGMM(pc, c_min, c_max, top = 2, plot = False):
    aic = []
    bic = []
    sil = []
    numberOfClusters = range (c_min, c_max)
    for n in numberOfClusters:
        model = GaussianMixture(n, covariance_type ='full', random_state = 0).fit(pc)
        clusters = model.predict(pc)
        bic.append(model.bic(pc))
        aic.append(model.aic(pc))
        sil.append(metrics.silhouette_score(pc,clusters ,metric='euclidean', sample_size=None, random_state=None))

    if plot:
        plt.plot(numberOfClusters, bic, label = 'BIC')
        plt.plot(numberOfClusters, aic, label = 'AIC')
        plt.legend()
        plt.title('BIC/AIC')
        plt.xlabel('n_components')
        plt.figure()
        plt.plot(numberOfClusters, sil, label = 'sil')
    
    bestBic = np.argsort(bic)[:top] + c_min
    bestAic = np.argsort(aic)[:top] + c_min
    bestSil = np.argsort(sil)[::-1][:top] + c_min
    return bestBic, bestAic, bestSil



def getUmap(dataset, pca_comp = 10):
    dataset = PCA(n_components=pca_comp).fit_transform(dataset)
    reducer = umap.UMAP(n_neighbors=50,
                          min_dist=0.1)
    embedding2d = reducer.fit_transform(dataset)
    return embedding2d


def plotBestPrediction(summaryDf, dataset, pca_comp = 10):
    df, truth = loadData(dataset)
    umap2D = getUmap(df, pca_comp = pca_comp)
    best = summaryDf[summaryDf['result'] == summaryDf['result'].min()].iloc[0].to_dict()
    _, clusters = run(best)
    plt.figure(figsize=(14, 5))
    plt.subplot(121)
    plt.title('Ground truth')
    plt.scatter(umap2D[:, 0], umap2D[:, 1], s = 4, c = truth.clusters)

    plt.subplot(122)
    plt.title(f"Best prediction rand index {-summaryDf['result'].min()}")
    plt.scatter(umap2D[:, 0], umap2D[:, 1], s = 4, c = clusters)
    
    
def run(params):
    df, truth = loadData(params['dataset'])
    
    # Preprocessing remove genes which don't appear in at least minCellsPerGene cells
    discreteDf = np.zeros(df.shape)
    discreteDf[np.where(df>0)] = 1
    genesToKeep = np.where(discreteDf.sum(axis = 0)>=params['minCellsPerGene'] )[0]
    df= df[df.columns[genesToKeep]]
    del discreteDf 
    
    # Remove genes which have a very low variance as they are expressed equally in all cells
    logDf = np.log1p(df)
    nonZeroMean = logDf.mean(axis = 0)
    nonZeroMean[nonZeroMean==0] = 1e-10
    dispersion =logDf.var(axis = 0)/nonZeroMean
    genesToKeep = np.where(dispersion>=params['minGeneDispersion'] )[0]
    df= df[df.columns[genesToKeep]]
    del logDf, nonZeroMean, dispersion

    if params['log']:
        df = np.log1p(df)
        
    # scaling
    if params['scaler'] == 'none':
        scaledDf = df.values
    if params['scaler'] == 'standardScaleGenes':
        scaledDf = StandardScaler().fit_transform(df)
    if params['scaler'] == 'standardScaleCells':
        scaledDf = StandardScaler().fit_transform(df.T).T
    if params['scaler'] == 'robustScaleGenes':
        scaledDf = RobustScaler().fit_transform(df)
    if params['scaler'] == 'robustScaleCells':
        scaledDf = RobustScaler().fit_transform(df.T).T
        
    # PCA reduction
    pc = PCA(n_components=params['pca_comp']).fit_transform(scaledDf)
    
    model = GaussianMixture(params['nb_clusters'] , covariance_type ='full', random_state = 0).fit(pc)
    clusters = model.predict(pc)
    score = adjusted_rand_score(truth.clusters.tolist(), clusters)
    params['randIndex'] = score
    return params, clusters