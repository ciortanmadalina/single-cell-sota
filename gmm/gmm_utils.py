import sys
sys.path.append("..") # this adds to path parent directory in order to import utils file

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
import itertools
from scipy.spatial.distance import cdist
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.neighbors import kneighbors_graph
import igraph as ig
import louvain
from sklearn.metrics.cluster import adjusted_rand_score
import umap
import os
from scipy import sparse, io

# Import utils
import utils
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
    best = summaryDf[summaryDf['_rand_index'] == summaryDf['_rand_index'].min()].iloc[0].to_dict()
    _, clusters = run(best)
    plt.figure(figsize=(14, 5))
    plt.subplot(121)
    plt.title('Ground truth')
    plt.scatter(umap2D[:, 0], umap2D[:, 1], s = 4, c = truth.clusters)

    plt.subplot(122)
    plt.title(f"Best prediction rand index {-summaryDf['_rand_index'].min()}")
    plt.scatter(umap2D[:, 0], umap2D[:, 1], s = 4, c = clusters)
    
    
def runGMM(params):
    df, truth = utils.loadData(params['dataset'])
    dfOrig = df.copy()
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
    ev = externalValidation(truth.clusters.tolist(), clusters)
    iv = internalValidation(dfOrig, clusters)

    params = {**params, **ev, **iv}
    return params, clusters



def cluster_knn_louvain(data, neighbors = 10):
    A = kneighbors_graph(data, 10, mode='connectivity', include_self=True)
    sources, targets = A.nonzero()
    weights = A[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=False)
    g.add_vertices(A.shape[0])  # this adds adjacency.shap[0] vertices
    g.add_edges(list(zip(sources, targets)))

    g.es['weight'] = weights
    weights = np.array(g.es["weight"]).astype(np.float64)
    partition_type = louvain.RBConfigurationVertexPartition
    partition_kwargs = {}
    partition_kwargs["weights"] = weights
    part = louvain.find_partition(g, partition_type, **partition_kwargs)
    groups = np.array(part.membership)
    return groups

def getUmap(data, ncomp = 2):
    reducer = umap.UMAP(n_neighbors=50,
                          min_dist=0.1, n_components=ncomp)
    embedding2d = reducer.fit_transform(data)
    return embedding2d

def runLouvain(params):
    df, truth = utils.loadData(params['dataset'])
    dfOrig = df.copy()
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
    data = PCA(n_components=params['pca_comp']).fit_transform(scaledDf)
    if params['doUmap']:
        data = getUmap(data, ncomp = params['umap_comp'])
    
    clusters = cluster_knn_louvain(data, neighbors = params['nb_neighbors'])
    ev = externalValidation(truth.clusters.tolist(), clusters)
    iv = internalValidation(dfOrig, clusters)

    params = {**params, **ev, **iv}
    return params, clusters


def externalValidation(truthClusters, predictedClusters):
    def purity_score(y_true, y_pred):
        # compute contingency matrix (also called confusion matrix)
        contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
        # return purity
        return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix) 
    scores = {}
    scores['_rand_index'] = adjusted_rand_score(truthClusters, predictedClusters)
    scores['_homogeneity_score'] = metrics.homogeneity_score(truthClusters, predictedClusters)
    scores['_purity_score'] = purity_score(truthClusters, predictedClusters)
    scores['_adjusted_mutual_info_score'] = metrics.adjusted_mutual_info_score(truthClusters, predictedClusters)
    scores['_fowlkes_mallows_score'] = metrics.fowlkes_mallows_score(truthClusters, predictedClusters)  
    return scores


def internalValidation(data, clusters):
    scores = {}
    """
    The score is bounded between -1 for incorrect clustering and +1 for highly dense clustering. 
    Scores around zero indicate overlapping clusters.
    The score is higher when clusters are dense and well separated, which relates to a standard concept of a cluster.
    """
    scores['_silhouette_score'] =metrics.silhouette_score(data,clusters ,metric='euclidean')
    """
    The score is higher when clusters are dense and well separated, which relates to a standard concept of a cluster.
    The score is fast to compute
    """
    scores['_calinski_harabaz_score'] = metrics.calinski_harabaz_score(data,clusters)
    """
    Zero is the lowest possible score. Values closer to zero indicate a better partition.
    The Davies-Boulding index is generally higher for convex clusters than other concepts of clusters, 
    such as density based clusters like those obtained from DBSCAN.
    """
    scores['_davies_bouldin_score'] = metrics.davies_bouldin_score(data,clusters)
    return scores


def plotCorrelation(resultsDf, name = None):
    scoreColumns = [c for c in resultsDf.columns if c.startswith('_')]
    score = resultsDf[scoreColumns]
    score = score.astype(float)
    if name is not None:
        plt.title(f"Correlation for {name}")
    sns.heatmap(score.corr(), annot=True)
