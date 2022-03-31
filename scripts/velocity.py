import numpy as np
import scvelo as scv
import pandas as pd
import loompy
import glob
import os
from scipy import io
import matplotlib.pylab as plt
import scanpy as sc
import logging


logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

logger.info('Reading loom file')
# loom_data = scv.read('/home/hsarkar/Projects/ThMYCN/tumor_merged.loom')
# gene_list = \
#         loom_data.var.index.str.split('_').str[1].str.upper().tolist()[:-1] + ['MYCN']

# loom_data.var.index = gene_list
# obs = loom_data.obs.copy()
# obs['aux_pre'] = obs.index.str.split(':').str[0]
# obs['cell_barcode'] = obs.index.str.split(':').str[1].str[:-1]
# obs['cell_id'] = obs.aux_pre + '_' + obs.cell_barcode + '-1'
# loom_data.obs.index = obs.cell_id.values
# logger.info('Created loom object')

loom_data = sc.read_h5ad(
    os.path.join(
        '/home/hsarkar/Projects',
        'neuroblastoma_analysis/results/scanpy',
        'loom_data.h5'
    )
)

cc_gene_file = var_gene_file = os.path.join(
    '/home/hsarkar/Projects/neuroblastoma_analysis/results/seurat_2',
    'cc_genes.csv'
) 

with open(cc_gene_file, 'r') as f:
    cell_cycle_genes = f.read().splitlines()

def run_velocity(out_dir):
    X1 = io.mmread(
    os.path.join(
        out_dir,
        'counts.mtx'
        )
    )
    logger.info('Construct adata...')
    adata = sc.AnnData(
            X = X1.transpose().tocsr()
    )

    gene_file = os.path.join(
        out_dir,
        'genes.csv'
    )
    logger.info('read genes...')
    with open(gene_file, 'r') as f:
        gene_names = f.read().splitlines()

    # var_gene_file = os.path.join(
    #     out_dir,
    #     'var_genes.csv'
    # )
    # logger.info('read genes...')
    # with open(var_gene_file, 'r') as f:
    #     var_gene_names = f.read().splitlines()

    logger.info('read metadata...')
    meta_data = pd.read_csv(
        os.path.join(
            out_dir,
            'meta_data.csv'
        ),
        index_col = 0
    )

    adata.obs = meta_data
    adata.var.index = gene_names

    pca = pd.read_csv(
        os.path.join(
            out_dir,
            'pca.csv'
        )
    )
    pca.index = adata.obs.index
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_umap'] = adata.obs[['UMAP_1', 'UMAP_2']].values 
    logger.info('merging velocity data...')
    adata_merged = scv.utils.merge(adata, loom_data)
    adata_merged.obs['clusters'] = adata_merged.obs.new_clusters

    logger.info('copying data for subsetting')
    adata_merged = adata_merged[:, 
                [x not in cell_cycle_genes for x in adata_merged.var.index]
    ].copy()
    adata_merged_backup = adata_merged.copy()


    # calculate velocity for the whole data
    # logger.info('Running velocity for the whole dataset')
    # scv.pp.filter_and_normalize(
    #     adata_merged, 
    #     min_shared_counts=20, 
    #     n_top_genes=2000
    # )
    # scv.tl.velocity(
    #     adata_merged, 
    #     mode='stochastic'
    # )

    # scv.tl.velocity_graph(
    #     adata_merged,
    #     n_jobs=30
    # )

    # # save the object
    # adata_merged.write_h5ad(
    #     os.path.join(
    #         out_dir,
    #         'neuroblast_5_samples_loom_merged.h5'
    #     )
    # )

    # take subset
    cur_celltypes = [
        'Immature ADRN',
        'Sympathoblasts 1',
        'Sympathoblasts 2',
        'Sympathoblasts 3',
        'Sympathoblasts 4',
        'Sympathoblasts 5',
        'Chromaffin cells'
    ]
    #'Sympathoblasts 6',
    
    logger.info('copying subset')
    adata_subset = adata_merged_backup[
        adata_merged_backup.obs['clusters'
    ].isin(cur_celltypes),].copy()
    logger.info('recompute umap')
    sc.pp.neighbors(adata_subset, n_neighbors = 30)
    sc.tl.umap(adata_subset)
    scv.pp.filter_and_normalize(
        adata_subset, 
        min_shared_counts=20, 
        n_top_genes=2000
    )
    # rerun velocity
    scv.tl.velocity(
        adata_subset, 
        mode='stochastic'
    )

    scv.tl.velocity_graph(
        adata_subset,
        n_jobs=30
    )

    adata_subset.write_h5ad(
        os.path.join(
            out_dir,
            'neuroblast_5_samples_subset_loom_merged.h5'
        )
    )

logger.info('cc_difference')
run_velocity('/home/hsarkar/Projects/neuroblastoma_analysis/results/seurat_2/cc_difference')

# logger.info('cc_remove')
# run_velocity('/home/hsarkar/Projects/neuroblastoma_analysis/results/seurat_2/cc_remove')