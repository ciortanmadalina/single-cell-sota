import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import pickle
import random
from tqdm import tqdm
import numpy as np
import statsmodels.api as sm

from sklearn.manifold import MDS
from mpl_toolkits import mplot3d
from scipy.spatial import distance
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering, SpectralClustering, KMeans, AffinityPropagation, DBSCAN, FeatureAgglomeration
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from IPython.display import clear_output, Image, display
from sklearn.datasets.samples_generator import make_blobs
import itertools
from scipy.spatial.distance import cdist
from sklearn.mixture import GaussianMixture
from mpl_toolkits.mplot3d import Axes3D
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
    
if printFunctionNames:
    print('bicAicAnalysis')
def bicAicAnalysis(X, numberOfClusters):
    print('Best score corresponds to the minimum values of AIC/BIC')
    bic = []
    aic = []

    for n in tqdm(numberOfClusters):
        model = GaussianMixture(n, covariance_type ='full', random_state = 0).fit(X)
        bic.append(model.bic(X))
        aic.append(model.aic(X))

    plt.plot(numberOfClusters, bic, label = 'BIC')
    plt.plot(numberOfClusters, aic, label = 'AIC')
    plt.legend()
    plt.title('BIC/AIC')
    plt.xlabel('n_components')