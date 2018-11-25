import numpy as np
import pandas as pd
import scanpy.api as sc
from scipy import sparse, io
from collections import Counter
import os.path
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib.pyplot as plt
plt.ion()
plt.show()
sc.settings.set_figure_params(dpi=80)
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = '../write/pbmc3k.h5ad'

def loadData(inputDataset):
    """
    Load input dataset
    """
    def loadAdata(path):
        adata = sc.read(f'{path}matrix.mtx', cache=True).T  # transpose the data
        adata.var_names = np.load(f'{path}/genes.npy')
        adata.obs_names = np.load(f'{path}/cells.npy')
        return adata
    
    if inputDataset == 'brainCIDR':
        path = '../input/brainCIDR/'
        if os.path.isfile(f'{path}matrix.mtx') == False:
            df = pd.read_csv(f"{path}brainTags.csv", index_col = 0)
            np.save(f'{path}cells.npy',  df.columns)
            np.save(f'{path}genes.npy', df.index.values)
            print(df.shape)
            df.head()
            m = sparse.csr_matrix(df.as_matrix())
            io.mmwrite(f'{path}matrix.mtx', m)
            del m, df
        adata = loadAdata(path)
        
        if os.path.isfile(f'{path}truth.pkl'):
            truth = pd.read_pickle(f'{path}truth.pkl')
        else:
            truth = pd.read_csv(f'{path}truth.csv', index_col = 0)
            truth.set_index('colnames.brainTags.', inplace=True)
            truth['truth'] = truth['truth'].astype('category')
            truth['clusters'] = truth['truth'].cat.codes
            truth.to_pickle(f'{path}truth.pkl')
        return adata, truth
    
    if inputDataset == 'pancreaticIsletCIDR':
        path = '../input/pancreaticIsletCIDR/'
        if os.path.isfile(f'{path}matrix.mtx') == False:
            df = pd.read_csv(f"{path}pancreaticIsletTags.csv", index_col = 0)
            np.save(f'{path}cells.npy',  df.columns)
            np.save(f'{path}genes.npy', df.index.values)
            print(df.shape)
            df.head()
            m = sparse.csr_matrix(df.as_matrix())
            io.mmwrite(f'{path}matrix.mtx', m)
            del m, df
        adata = loadAdata(path)
        
        if os.path.isfile(f'{path}truth.pkl'):
            truth = pd.read_pickle(f'{path}truth.pkl')
        else:
            truth = pd.read_csv(f'{path}truth.csv', index_col = 0)
            truth.set_index('colnames.pancreaticIsletTags.', inplace=True)
            truth['truth'] = truth['truth'].astype('category')
            truth['clusters'] = truth['truth'].cat.codes
            truth.to_pickle(f'{path}truth.pkl')
        return adata, truth
    
    
    if inputDataset == 'deng':
        path = '../input/deng/'
        if os.path.isfile(f'{path}matrix.mtx') == False:
            df = pd.read_csv(f"{path}deng.csv", index_col = 0)
            np.save(f'{path}cells.npy',  df.columns)
            np.save(f'{path}genes.npy', df.index.values)
            print(df.shape)
            df.head()
            m = sparse.csr_matrix(df.as_matrix())
            io.mmwrite(f'{path}matrix.mtx', m)
            del m, df
        adata = loadAdata(path)
        
        if os.path.isfile(f'{path}truth.pkl'):
            truth = pd.read_pickle(f'{path}truth.pkl')
        else:
            truth = pd.read_csv(f'{path}truth.csv', index_col = 0)
            truth['cells'] = np.load(f'{path}cells.npy')
            truth.set_index('cells', inplace=True)
            truth['truth'] = truth['x'].astype('category')
            truth['clusters'] = truth['truth'].cat.codes
            truth.drop('x', axis = 1, inplace = True)
            truth.to_pickle(f'{path}truth.pkl')
        return adata, truth
    
def preprocess(adata, min_genes, min_cells, teta_total_features, normalize_per_cell,
               filter_min_mean, filter_max_mean, filter_min_disp, regress_out, scale ):
    """
    Preprocessing phase
    """

    adata.var_names_make_unique()
    
    if min_genes is not None: 
        sc.pp.filter_cells(adata, min_genes=min_genes)  
        
    if min_cells is not None: 
        sc.pp.filter_genes(adata, min_cells=min_cells)
        
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    sc.pl.scatter(adata, x='total_counts', y='total_features_by_counts')
    if teta_total_features is not None: 
        adata = adata[adata.obs['total_features_by_counts'] < teta_total_features, :]
        
    adata.raw = sc.pp.log1p(adata, copy=True)
    
    if normalize_per_cell :
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    
    if filter_min_mean is not None:
        filter_result = sc.pp.filter_genes_dispersion(
                        adata.X, min_mean=filter_min_mean, 
                        max_mean=filter_max_mean, 
                        min_disp=filter_min_disp)
        print(f'Filtering {Counter(filter_result.gene_subset)}')
        sc.pl.filter_genes_dispersion(filter_result)
        
        adata = adata[:, filter_result.gene_subset]
        print(f'Keeping {len(adata.var_names)} genes')
    else:
        print('No dispertion gene filter')
    sc.pp.log1p(adata)
    
    if regress_out is not None:
        sc.pp.regress_out(adata, [regress_out])
    if scale is not None:
        sc.pp.scale(adata, max_value=scale)
    return adata


def clustering(adata, n_neighbors, n_pcs, plot_pca):
    """
    Cluster data
    """
    sc.tl.pca(adata, svd_solver='arpack')
    if plot_pca:
        sc.pl.pca(adata)
        sc.pl.pca_variance_ratio(adata, log=True)
    sc.pp.neighbors(adata, n_neighbors= n_neighbors, n_pcs= n_pcs)
    sc.tl.umap(adata)
    sc.tl.louvain(adata)
    clusters = adata.obs['louvain'].get_values()

    return adata.obs_names, clusters


def evaluate(truth, cells, clusters, umapCoord):
    score = adjusted_rand_score(
                truth[truth.index.isin(cells)].clusters.tolist(), 
                clusters)
    if umapCoord is not None:
        plt.figure(figsize = (12, 5))
        
        plt.subplot(121)
        plt.title('Ground truth')
        plt.scatter(umapCoord[:,0], 
                    umapCoord[:,1], 
                    c = truth[truth.index.isin(cells)].clusters.tolist())
        
        
        plt.subplot(122)
        plt.title('Clustering')
        plt.scatter(umapCoord[:,0], 
                    umapCoord[:,1], 
                    c = clusters)
    return score

def filterDictionary(params, prefix = 'preprocess_'):
    return {k[len(prefix):] : v for k, v in params.items() if k.startswith(prefix)}

def run(params):
    adata, truth = loadData(params['load_inputDataset'])
    print(f"Loading dataset {params['load_inputDataset']} with {adata.var_names.shape[0]} genes and {adata.obs_names.shape[0]} cells")

    preprocess(adata,**filterDictionary(params, prefix = 'preprocess_'))
    cells, clusters = clustering(adata, **filterDictionary(params, prefix = 'cluster_'))
    evalParams = {}
    evalParams['cells'] = cells
    evalParams['clusters'] = clusters
    evalParams['umapCoord'] = adata.obsm['X_umap'] if params['evaluate_plot_results'] else None
    randIndex = evaluate(truth, **evalParams)
    params['randIndex'] = randIndex
    print (f"Rand_index {randIndex}")
    return params