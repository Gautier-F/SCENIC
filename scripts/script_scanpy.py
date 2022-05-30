import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE

import seaborn as sns
import matplotlib.pyplot as plt


# f_loom_path_unfilt = "loomfiles/dat_unfilt.loom"
f_loom_path_scenic = snakemake.output[0]
# f_anndata_path = "loomfiles/anndata.h5ad"


sc.settings.njobs=20

f_mtx_dir = snakemake.input[0]
adata = sc.read_10x_mtx(f_mtx_dir, var_names='gene_symbols', cache= False)

# row_attrs = { 
    # "Gene": np.array(adata.var.index) ,
# }
# col_attrs = { 
    # "CellID":  np.array(adata.obs.index) ,
    # "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    # "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
# }

# lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs )

#####QC
sc.pp.filter_cells(adata, min_genes=0)
mito_genes=adata.var_names.str.startswith('MT-')
adata.obs['percent_mito']=np.sum(adata[:,mito_genes].X, axis=1).A1/np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)

sns.distplot( adata.obs['n_genes'], ax=ax1, norm_hist=True, bins=100)
sns.distplot( adata.obs['n_counts'], ax=ax2, norm_hist=True, bins=100)
sns.distplot( adata.obs['percent_mito'], ax=ax3, norm_hist=True, bins=100)

ax1.title.set_text('Number of genes expressed per cell')
ax2.title.set_text('Counts per cell')
ax3.title.set_text('Mitochondrial read fraction per cell')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

fig.tight_layout()

fig.savefig(snakemake.output['pdf1'], dpi=600, bbox_inches='tight')

### filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells = 3)
adata = adata[adata.obs['n_genes'] <6000, :]
adata = adata[adata.obs['percent_mito'] < 0.2, :]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)

adata.obs['n_genes']

sns.distplot( adata.obs['n_genes'], ax=ax1, norm_hist=True, bins=100)
sns.distplot( adata.obs['n_counts'], ax=ax2, norm_hist=True, bins=100)
sns.distplot( adata.obs['percent_mito'], ax=ax3, norm_hist=True, bins=100)

ax1.title.set_text('Number of genes expressed per cell')
ax2.title.set_text('Counts per cell')
ax3.title.set_text('Mitochondrial read fraction per cell')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

fig.tight_layout()

fig.savefig(snakemake.output['pdf2'], dpi=600, bbox_inches='tight')

# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( snakemake.output['loom'], adata.X.transpose(), row_attrs, col_attrs)

nGenesDetectedPerCell = np.sum(adata.X>0, axis=1)
nGenesDetectedPerCell = pd.DataFrame(nGenesDetectedPerCell)
percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])

fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
for i,x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
    ax.text(x=x, y=ax.get_ylim()[1], s=f"{int(x)} ({percentiles.index.values[i]*100}%)", color='red', rotation=30, size='x-small',rotation_mode='anchor' )
ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')
fig.tight_layout()
fig.savefig(snakemake.output['pdf3'], dpi=600, bbox_inches='tight')

# adata.write(f_anndata_path)

























